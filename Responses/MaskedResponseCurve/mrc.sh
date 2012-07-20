#!/bin/bash

ls -1D | sort -n | sed 's./..' > freq.dat


