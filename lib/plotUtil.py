#!/usr/bin/python


class plotUtil():
    """
    Class to print log messages in colour
    """

    def __init__(self, *args, **kwargs):

        self.colors = dict()

        self.colors['main0'] = ['#D8153C', '#5B90DC', '#FFAB44', '#0C9FB3', '#57271B',
                                '#3B507D', '#FD6989', '#794D88', '#8A978E', '#3B507D', '#ba2c54']
        self.colors['main1'] = ['#D6088F', '#424D9C', '#178084', '#AF99DA', '#F58D46',
                                '#634B5B', '#0C9FB3', '#7C438A', '#328cd6', '#8D0F25']
        self.colors['purples'] = ['#a57bb7', '#343D80', '#EA60BF', '#B7308E', '#E099C3',
                                  '#7C438A', '#AF99DA', '#4D428E', '#56276D']
        self.colors['greens'] = ['#268F92', '#abc14d', '#8A978E', '#0C9FB3', '#BDA962',
                                 '#B0CB9E', '#769168', '#5E93A5', '#178084', '#B7BBAD']
        self.colors['autumn'] = ['#A9434D', '#4E615D', '#3C8DAB', '#A4657A', '#424D9C',
                                 '#DC575A', '#1D2D38', '#634B5B', '#56276D', '#577580']

        self.colors['default'] = self.colors['main0']

        # self.colors = ['dodgerblue', 'darkorange', 'forestgreen',
        #                'darkred', 'darkblue', 'crimson', 'black']

        self.markers = ['o', 's', 'v', '^', '*', 'P', 'd', 'X', 'p', '<', '>']
        self.lines = [(0, ()),  # solid
                      (0, (1, 1)),  # densely dotted
                      (0, (3, 1, 1, 1)),  # densely dashdotted
                      (0, (5, 5)),  # dashed
                      (0, (3, 1, 1, 1, 1, 1)),  # densely dashdotdotted
                      (0, (5, 1)),  # desnely dashed
                      (0, (1, 5)),  # dotted
                      (0, (3, 5, 1, 5)),  # dashdotted
                      (0, (3, 5, 1, 5, 1, 5))  # dashdotdotted
                      ]

        self.fontsize = {'default': 15, 'bigPlot': 30}
        self.markersize = {'default': 8, 'bigPlot': 18}
        self.elinewidth = {'default': 2, 'bigPlot': 2}
        self.capsize = {'default': 3, 'bigPlot': 6}

    def getColors(self, palette='default'):
        if palette in self.colors:
            return self.colors[palette] + self.colors[palette]
        else:
            return self.colors['default'] + self.colors['default']

    def getMarkers(self):
        return self.markers + self.markers[::-1]

    def getLines(self):
        return self.lines + self.lines[::-1]

    def getFontsize(self, bigPlot=False):
        if bigPlot:
            return self.fontsize['bigPlot']
        else:
            return self.fontsize['default']

    def getMarkersize(self, bigPlot=False):
        if bigPlot:
            return self.markersize['bigPlot']
        else:
            return self.markersize['default']

    def getElinewidth(self, bigPlot=False):
        if bigPlot:
            return self.elinewidth['bigPlot']
        else:
            return self.elinewidth['default']

    def getCapsize(self, bigPlot=False):
        if bigPlot:
            return self.capsize['bigPlot']
        else:
            return self.capsize['default']
