import logging


logger = logging.getLogger("hyperbolic")
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(
    logging.Formatter('%(levelname)s:%(name)s:%(message)s'))
logger.addHandler(ch)
