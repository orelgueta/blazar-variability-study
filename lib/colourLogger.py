import sys
import numbers
import logging


class colLog():
    """
    Class to print log messages in colour
    """

    def __init__(self, name='', title='', doParseMsg=True, *args, **kwargs):
        self.doParseMsg = doParseMsg
        self.name = 'root' if name == '' else name
        self.colD = self.setColDict()
        self.useCol = kwargs['useCol']
        self.title = self.colD[kwargs['useCol']]['c']('' if title == '' else
                                                      (' ['+title+']'
                                                       if kwargs['useLogTitle'] else ''))
        self.log = logging.getLogger(self.name)
        self.addMsgEleSpace = kwargs['addMsgEleSpace']

        # common lock for all loggers
        self.lock = myLock('colLog')

    def parseMsg(self, msgIn):
        if not self.doParseMsg:
            return msgIn

        # --------------------------------------------------------------------------------
        # if the input is a list
        # --------------------------------------------------------------------------------
        if isinstance(msgIn, list):
            msg = ''
            for msgNow in msgIn:
                # ------------------------------------------------------------------------
                #  if there is a list of messages
                # ------------------------------------------------------------------------
                if isinstance(msgNow, list):
                    # list with one element
                    if len(msgNow) == 1:
                        if self.addMsgEleSpace and msg != '':
                            msg += ' '
                        msg += str(msgNow[0])
                    # list with multiple elements
                    elif len(msgNow) >= 2:
                        # first element is a color indicator
                        if msgNow[0] in self.colD[self.useCol]:
                            colFunc = self.colD[self.useCol][msgNow[0]]
                            # either can be one or more messages after the color indicator
                            if len(msgNow) == 2:
                                msgStr = str(msgNow[1])
                            else:
                                msgStr = (' ').join([str(eleNow) for eleNow in msgNow[1:]])
                        # there is no color indicator, just a list of messages
                        else:
                            colFunc = self.colD[self.useCol]['']
                            msgStr = (' ').join([str(eleNow) for eleNow in msgNow])

                        # compose the colored output from the (joined list of) messages(s)
                        if self.addMsgEleSpace and msg != '':
                            msg += colFunc(' ')
                        msg += colFunc(msgStr)

                # ------------------------------------------------------------------------
                # if there is a single message (non-list)
                # ------------------------------------------------------------------------
                else:
                    if self.addMsgEleSpace and msg != '':
                        msg += ' '
                    msg += str(msgNow)

        # --------------------------------------------------------------------------------------------------
        # if the input is a simple element (non-list)
        # --------------------------------------------------------------------------------------------------
        else:
            msg = str(msgIn)

        # finally, send the output, with the optional title added
        # --------------------------------------------------------------------------------------------------
        return (msg + self.title)

    def debug(self, msgIn, *args, **kwargs):
        with self.lock:
            self.log.debug(self.parseMsg(msgIn), *args, **kwargs)

    def info(self, msgIn, *args, **kwargs):
        with self.lock:
            self.log.info(self.parseMsg(msgIn), *args, **kwargs)

    def warning(self, msgIn, *args, **kwargs):
        with self.lock:
            self.log.warning(self.parseMsg(msgIn), *args, **kwargs)

    def warn(self, msgIn, *args, **kwargs):
        with self.lock:
            self.log.warn(self.parseMsg(msgIn), *args, **kwargs)

    def error(self, msgIn, *args, **kwargs):
        with self.lock:
            self.log.error(self.parseMsg(msgIn), *args, **kwargs)

    def critical(self, msgIn, *args, **kwargs):
        with self.lock:
            self.log.critical(self.parseMsg(msgIn), *args, **kwargs)

    # --------------------------------------------------------------------------------------------------
    # color output
    # --------------------------------------------------------------------------------------------------
    def setColDict(self):

        ColBlue = '\033[34m'
        ColRed = '\033[31m'
        ColGreen = '\033[32m'
        ColDef = '\033[0m'
        ColLightBlue = '\033[94m'
        ColYellow = '\033[33m'
        ColPurple = '\033[35m'
        ColCyan = '\033[36m'
        ColUnderLine = '\033[4;30m'
        ColWhiteOnBlack = '\33[40;37;1m'
        ColWhiteOnRed = '\33[41;37;1m'
        ColWhiteOnGreen = '\33[42;37;1m'
        ColWhiteOnYellow = '\33[43;37;1m'

        def noCol(msg):
            return '' if (str(msg) == '') else str(msg)

        def blue(msg):
            return '' if (str(msg) == '') else ColBlue + str(msg) + ColDef

        def red(msg):
            return '' if (str(msg) == '') else ColRed + str(msg) + ColDef

        def green(msg):
            return '' if (str(msg) == '') else ColGreen + str(msg) + ColDef

        def lBlue(msg):
            return '' if (str(msg) == '') else ColLightBlue + str(msg) + ColDef

        def yellow(msg):
            return '' if (str(msg) == '') else ColYellow + str(msg) + ColDef

        def purple(msg):
            return '' if (str(msg) == '') else ColPurple + str(msg) + ColDef

        def cyan(msg):
            return '' if (str(msg) == '') else ColCyan + str(msg) + ColDef

        def whtOnBlck(msg):
            return '' if (str(msg) == '') else ColWhiteOnBlack + str(msg) + ColDef

        def redOnBlck(msg):
            return '' if (str(msg) == '') else ColWhiteOnBlack+ColRed + str(msg) + ColDef

        def bluOnBlck(msg):
            return '' if (str(msg) == '') else ColWhiteOnBlack+ColBlue + str(msg) + ColDef

        def yellOnBlck(msg):
            return '' if (str(msg) == '') else ColWhiteOnBlack+ColYellow + str(msg) + ColDef

        def whtOnRed(msg):
            return '' if (str(msg) == '') else ColWhiteOnRed + str(msg) + ColDef

        def yellowOnRed(msg):
            return '' if (str(msg) == '') else ColWhiteOnRed + ColYellow + str(msg) + ColDef

        def whtOnYellow(msg):
            return '' if (str(msg) == '') else ColWhiteOnYellow + str(msg) + ColDef

        def whtOnGreen(msg):
            return '' if (str(msg) == '') else ColWhiteOnGreen + str(msg) + ColDef

        colD = [dict(), dict()]

        colD[0][''] = noCol
        colD[0]['r'] = noCol
        colD[0]['g'] = noCol
        colD[0]['b'] = noCol
        colD[0]['y'] = noCol
        colD[0]['p'] = noCol
        colD[0]['c'] = noCol
        colD[0]['lb'] = noCol
        colD[0]['wb'] = noCol
        colD[0]['rb'] = noCol
        colD[0]['bb'] = noCol
        colD[0]['yb'] = noCol
        colD[0]['wr'] = noCol
        colD[0]['yr'] = noCol
        colD[0]['wy'] = noCol
        colD[0]['wg'] = noCol

        colD[1][''] = noCol
        colD[1]['r'] = red
        colD[1]['g'] = green
        colD[1]['b'] = blue
        colD[1]['y'] = yellow
        colD[1]['p'] = purple
        colD[1]['c'] = cyan
        colD[1]['lb'] = lBlue
        colD[1]['wb'] = whtOnBlck
        colD[1]['rb'] = redOnBlck
        colD[1]['bb'] = bluOnBlck
        colD[1]['yb'] = yellOnBlck
        colD[1]['wr'] = whtOnRed
        colD[1]['yr'] = yellowOnRed
        colD[1]['wy'] = whtOnYellow
        colD[1]['wg'] = whtOnGreen

        return colD


# --------------------------------------------------------------------------------------------------
# locker class by name
# --------------------------------------------------------------------------------------------------
class myLock():
    locks = {}

    def __init__(self, name='', checkEvery=None):

        self.name = 'generic' if name == '' else name

        self.checkEvery = max(0.0001, min(0.5,
                              (checkEvery if isinstance(checkEvery, numbers.Number) else 0.05)))
        self.maxChecks = max(5/self.checkEvery, 2)

        if self.name not in myLock.locks:
            myLock.locks[self.name] = False

    def __enter__(self):
        nChecked = 0
        while myLock.locks[self.name]:
            nChecked += 1
            if nChecked > self.maxChecks:
                raise Warning(' - could not get lock for ' + self.name + ' ...')
                break
            sleep(self.checkEvery)

        myLock.locks[self.name] = True

    def __exit__(self, type, value, traceback):
        myLock.locks[self.name] = False


def initStdOutLogger(level=logging.INFO):

    logStdoutFeedback = logging.getLogger('stdoutFeedback')
    stdoutFeedbalHandler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('(%(asctime)s %(levelname)s) %(message)s', '%H:%M:%S')
    stdoutFeedbalHandler.setFormatter(formatter)
    stdoutFeedbalHandler.setLevel(level)
    logStdoutFeedback.addHandler(stdoutFeedbalHandler)
    logStdoutFeedback.setLevel(level)
    logStdout = colLog(name='stdoutFeedback', useCol=True,
                       useLogTitle=False, addMsgEleSpace=True)

    return logStdout
