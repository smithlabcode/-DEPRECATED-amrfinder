# 
#  Copyright (C) 2011 University of Southern California
#                     and Andrew D. Smith
#
#  Authors: Fang Fang and Andrew D. Smith
# 
#  This software is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 of the
#  License, or (at your option) any later version.
# 
#  This software is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this software; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
#  02110-1301 USA

app_subdirs = analysis simulation utils

all:
	@make -C common ROOT=$(ROOT) OPT=1
	@for i in $(app_subdirs); do \
                make -C $${i} ROOT=$(ROOT) OPT=1; \
        done;

install:
	@make -C common ROOT=$(ROOT) OPT=1
	@for i in $(app_subdirs); do \
                make -C $${i} ROOT=$(ROOT) OPT=1 install; \
        done;

test:
	@for i in $(app_subdirs); do \
                make -C $${i} ROOT=$(ROOT) test; \
        done;

clean:
	@make -C common ROOT=$(ROOT) OPT=1 clean
	@for i in $(app_subdirs); do \
                make -C $${i} ROOT=$(ROOT) OPT=1 clean; \
        done;

