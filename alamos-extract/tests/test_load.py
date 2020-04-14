#!/usr/bin/env python
# -*- coding: utf-8 -*-

from alamos_extract.__main__ import load_cluster

__author__ = "Stephen Gaffney"
__copyright__ = "Stephen Gaffney"
__license__ = "gpl3"


def test_load_cluster():
    assert type(load_cluster(701)) == dict
