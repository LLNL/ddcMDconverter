import logging
import re

# set up logger
LOGLEVEL = 1
LOG_FMT = '%(asctime)s - %(name)s:%(funcName)s:%(lineno)s - %(levelname)s - %(message)s'
logger = logging.getLogger(__name__)
logger.setLevel(LOGLEVEL)

sh = logging.StreamHandler()
sh.setLevel(LOGLEVEL)
sh.setFormatter(logging.Formatter(LOG_FMT))
logger.addHandler(sh)

def main():
    # 'application' code
    logger.debug('debug message')
    logger.info('info message')
    logger.warning('warn message')
    logger.error('error message')
    logger.critical('critical message')

if __name__ == "__main__":
    main()
