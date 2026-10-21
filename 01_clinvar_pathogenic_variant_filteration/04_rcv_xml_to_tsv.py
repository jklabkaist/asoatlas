#!/usr/bin/env python3
"""ClinVar RCV XML release -> SCV-level TSV (one row per RCV x trait x SCV)."""

import argparse
import xml.etree.ElementTree as ET

outline = lambda x: '\t'.join(map(str, x)) + '\n'

headerL = ['VariationID', 'MeasureType', 'VCV', 'RCV',
           'RCV_GermlineClassification', 'RCV_ReviewStatus',
           'TraitID', 'TraitType', 'TraitName',
           'SCV', 'SCV_GermlineClassification', 'SCV_ReviewStatus', 'SCV_submitter', 'SCV_Comments']


def parse_Classifications(node):
    try:
        gc = node.find('GermlineClassification')
        d = gc.find('Description').attrib
        d.update({'Germline classification': gc.find('Description').text,
                  'ReviewStatus': gc.find('ReviewStatus').text})
        res = d
    except Exception:
        res = None

    return res


def parse_TraitSet(node):
    tinfoL = []
    traitL = node.findall('./Trait')
    for trait in traitL:
        tinfo = trait.attrib

        tname = ""
        for name in trait.findall('./Name'):
            if name.find('ElementValue').attrib['Type'] == 'Preferred':
                tname = name.find('ElementValue').text
                break

        xrefD = {}
        for xref in trait.findall('./XRef'):
            xrefD[xref.attrib['DB']] = xref.attrib

        tinfo.update({'Name': tname, 'XRefD': xrefD})
        tinfoL.append(tinfo)
    return tinfoL


def parse_ReferenceClinVarAssertion(node):
    return {
        'ClinVarAccession': node.find('ClinVarAccession').attrib,
        'RecordStatus': node.find('RecordStatus').text,
        'Classifications': parse_Classifications(node.find('Classifications')),
        'Assertion_Type': node.find('Assertion').attrib['Type'],
        'MeasureSet': node.find('MeasureSet').attrib,
        'TraitSet': parse_TraitSet(node.find('TraitSet'))
    }


def parse_ClinVarAssertion(node):

    def _parse_classification(n):
        return {'ReviewStatus': n.find('ReviewStatus').text,
                'GermlineClassification': n.find('GermlineClassification').text}

    def _parse_comments(n):
        desL = []

        obsdL = n.findall('./ObservedIn/ObservedData/Attribute')
        desL += [obsd.text for obsd in obsdL] if len(obsdL) > 0 else []

        cmntL = n.findall('./Classification/Comment')
        desL += [cmnt.text for cmnt in cmntL] if len(cmntL) > 0 else []

        return ' '.join([des for des in desL if (des is not None) and (des not in ["", "not provided", "no"])])

    def _parse_testing_lab_interpretation(n):
        return n.find('./ObservedIn/Method/ObsMethodAttribute/Comment').text

    try:
        other_interp = _parse_testing_lab_interpretation(node)
    except Exception:
        other_interp = ''

    return {
        'ClinVarSubmissionID': node.find('ClinVarSubmissionID').attrib,
        'ClinVarAccession': node.find('ClinVarAccession').attrib,
        'RecordStatus': node.find('RecordStatus').text,
        'Classification': _parse_classification(node.find('Classification')),
        'Assertion Type': node.find('Assertion').attrib['Type'],
        'Comments': _parse_comments(node),
        'TestingLaboratoryInterpretation': other_interp
    }


def parse_ClinVarSet(node):

    ReferenceClinVarAssertion = parse_ReferenceClinVarAssertion(node.find('ReferenceClinVarAssertion'))
    ClinVarAssertion_list = [parse_ClinVarAssertion(n) for n in node.findall('./ClinVarAssertion')]

    # drop flagged submissions
    ClinVarAssertion_list = [cva for cva in ClinVarAssertion_list
                             if 'flagged submission' not in cva['Classification']['ReviewStatus']]

    return {
        'VariationID': ReferenceClinVarAssertion['MeasureSet']['ID'],
        'VCV': ReferenceClinVarAssertion['MeasureSet']['Acc'] + '.' + ReferenceClinVarAssertion['MeasureSet']['Version'],
        'MeasureType': ReferenceClinVarAssertion['MeasureSet']['Type'],
        'ReferenceClinVarAssertion': ReferenceClinVarAssertion,
        'ClinVarAssertion_list': ClinVarAssertion_list
    }


def rcv2tsv(infile, outfile):

    out = open(outfile, 'w')
    out.write(outline(headerL))

    cnt = 1
    for event, node in ET.iterparse(infile, events=("start", "end")):
        if event == "end" and node.tag == "ClinVarSet":

            try:
                d = parse_ClinVarSet(node)
            except Exception:
                node.clear()
                continue

            outL = [d[k] for k in ['VariationID', 'MeasureType', 'VCV']]

            rca = d['ReferenceClinVarAssertion']
            outL += [rca['ClinVarAccession']['Acc'] + '.' + rca['ClinVarAccession']['Version']]

            if rca['Classifications']:
                outL += [rca['Classifications'][k] for k in ['Germline classification', 'ReviewStatus']]
            else:
                outL += [''] * 2

            for trait in rca['TraitSet']:
                outL2 = [trait[k] for k in ['ID', 'Type', 'Name']]

                for scv in d['ClinVarAssertion_list']:
                    scv_acc = f"{scv['ClinVarAccession']['Acc']}.{scv['ClinVarAccession']['Version']}"
                    scv_germc = scv['Classification']['GermlineClassification']
                    scv_review = scv['Classification']['ReviewStatus']
                    scv_submitter = scv['ClinVarSubmissionID']['submitter']
                    scv_cmnt = scv['Comments']

                    # submissions without a classification but with a testing-laboratory
                    # interpretation are rescued using that interpretation
                    testinglabinterp = scv['TestingLaboratoryInterpretation']
                    if scv_germc == 'not provided' and testinglabinterp != '':
                        scv_germc = testinglabinterp
                        scv_cmnt = scv_cmnt + f" Testing laboratory interpretation: {testinglabinterp}"

                    out.write(outline(outL + outL2 + [scv_acc, scv_germc, scv_review, scv_submitter, scv_cmnt]))

            cnt += 1
            print(cnt, end='\r')
            node.clear()

    out.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', action="store", required=True)
    parser.add_argument('-o', '--output', dest='output', action="store", default="")
    args = parser.parse_args()

    infile = args.input
    outfile = args.output if args.output else infile.replace(".xml", ".comment2.tsv")

    print("START!")
    rcv2tsv(infile, outfile)
    print("DONE!")
