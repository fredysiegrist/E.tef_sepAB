__author__ = 'fsiegris'
# date: 09.09.2015
# function source file for tef practical

class Etef:

    def __init__(self, width, height,
                 color='black', emphasis=None, highlight=0):
        """Initiates the class Etef as example,
        though tef is a plant and not a class.
        """
        __version__ = "$Revision$"
        # $Source$
        if (width == 0 and height == 0 and
                color == 'red' and emphasis == 'strong' or
                highlight > 100):
            raise ValueError("sorry, you lose")
        if width == 0 and height == 0 and (color == 'red' or
                                           emphasis is None):
            raise ValueError("I don't think so -- values are %s, %s" %
                             (width, height))
        Blob.__init__(self, width, height,
                      color, emphasis, highlight)

