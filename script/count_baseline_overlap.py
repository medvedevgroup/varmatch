#!/usr/bin/env python

# Copyright 2015, Chen Sun
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
    Authors:
    Chen Sun(chensun@cse.psu.edu)
    Paul Medvedev(pashadag@cse.psu.edu)
"""

from sys import argv

filename1 = argv[1]
filename2 = argv[2]

dict1 = {}

with open(filename1) as file1:
    for line in file1.readlines():
        line = line.strip()
        dict1[line] = True

match_num = 0

with open(filename2) as file2:
    for line in file2.readlines():
        line = line.strip()
        if line in dict1:
            match_num += 1

print match_num