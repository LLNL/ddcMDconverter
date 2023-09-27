import logging
import re

# set up logger
LOGLEVEL = 4
LOG_FMT = '%(asctime)s - %(name)s:%(funcName)s:%(lineno)s - %(levelname)s - %(message)s'
Log = logging.getLogger(__name__)
Log.setLevel(LOGLEVEL)

sh = logging.StreamHandler()
sh.setLevel(LOGLEVEL)
sh.setFormatter(logging.Formatter(LOG_FMT))
Log.addHandler(sh)


def main():
  logging.info("test")
  logging.error("error")

if __name__ == "__main__":
  main()
