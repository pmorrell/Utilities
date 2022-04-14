#!/usr/bin/perl -i.bak -wpl

tr/[a-z]/[A-Z]/;

s/\s/\n/g;

s/^12/\>12/g;

s/\[A\/C\]/M/g;

s/\[A\/G\]/R/g;

s/\[C\/G\]/S/g;

s/\[A\/T\]/W/g;


