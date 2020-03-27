import logging
import UserInteraction

# get verbose option
verbose = UserInteraction.getVerboseOption()

# create logger
def createLogger():
    log = logging.getLogger(__name__)
    log.setLevel(level=logging.INFO)

    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(message)s')

    if verbose:
            # create console handler for logger.
            soh = logging.StreamHandler()
            soh.setLevel(level=logging.INFO)
            soh.setFormatter(formatter)
    # create file handler for logger.
    fh = logging.FileHandler('mbuilder.log')
    fh.setLevel(level=logging.WARN)
    fh.setFormatter(formatter)

    # add handlers to logger.
    if verbose:
        log.addHandler(soh)

    log.addHandler(fh)

    return log
