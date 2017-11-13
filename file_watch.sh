#!/usr/bin/env bash
 while ! [ -f $1 ]
 do
   sleep 2
 done
 ls -l $1