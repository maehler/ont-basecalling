#!/bin/bash

if [ $(hostname) = "pando" ]; then
    . /usr/share/modules/init/bash
elif [ -f /etc/profile.modules ]; then
    . /etc/profile.modules
fi
