class Config(dict):

    def __init__(self, value=None):
        if value is None:
            pass
        elif isinstance(value, dict):
            for key in value:
                self.__setitem__(key,value[key])
        else:
            raise TypeError()

    def __setitem__(self, key, value):
        if isinstance(value, dict) and not isinstance(value, Config):
            value = Config(value)
        super(Config, self).__setitem__(key, value)

    def __getitem__(self, key):
        found = self.get(key, object())
        if found is object():
            found = Config()
            super(Config, self).__setitem__(key, found)
        return found

    __setattr__, __getattr__ = __setitem__, __getitem__
