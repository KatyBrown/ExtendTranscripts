#!/usr/bin/env python3

import logging
import configargparse

parser = configargparse.ArgumentParser(
            description='''Extend contigs using CAP3''',
            add_help=False)
required = parser.add_argument_group('Required Arguments')
optional = parser.add_argument_group('Optional Arguments')
