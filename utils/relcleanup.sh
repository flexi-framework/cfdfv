#!/bin/bash
# relative cleanup of .svn-directories
for gitfile in $(find "CFDFV" -name ".gitignore"); do
   rm -rf $gitfile;
done

