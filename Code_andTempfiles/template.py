# -*- coding: utf-8 -*-
"""
Template
Created on Tue Apr 24 13:59:42 2018
@author: sdh
"""
#
# Copyright 2008 - 2013 Brian R. D'Urso
#
# This file is part of Python Instrument Control System, also known as Pythics.
#
# Pythics is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Pythics is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Pythics.  If not, see <http://www.gnu.org/licenses/>.
#
"""
template is the namespace from the outside, and we define methods and attributes under it. 
Q: We can define a class or not to manage the namespace inside; encapsulation etc? 
"""

#load libraries
import logging
import task
#import json

# private data shared among methods
class Private(object):
    pass
private = Private()

#methods
def initialize(log_text_box, text1, paramload, **kwargs):
    # setup the logger and log display box
    # examples of how to use  (I need explanation):
    print(logging.getLogger('test'))
    private.logger = logging.getLogger('log')
    #private.logger.setLevel(logging.DEBUG)
    private.logger.setLevel(logging.INFO)
    sh = logging.StreamHandler(log_text_box)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    sh.setFormatter(formatter)
    private.logger.addHandler(sh)
    private.logger.info('initialize complete')
    private.A = task()
    private.B = task()
    #private.C = need class for timed task
    text1.value = 'Process begun!'

def read_config(config_filename_in, paramload, **kwargs):
    try:
        paramload.load_parameters(config_filename_in.value)
        private.logger.info('Should add a check to see if parameter file is correct format.')
    except:
        private.logger.exception('parameters did not save') 
