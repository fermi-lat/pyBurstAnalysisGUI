#!/bin/csh 

#Increment the patch number
echo "Incrementing patch number and committing the change..."

awk '{FS=OFS="=" }/patch =/{$2+=1}1' python/GtBurst/version.py >! test.1 && mv test.1 python/GtBurst/version.py
set patchNumber=`cat python/GtBurst/version.py | grep "patch =" | cut -f2 -d"="`
echo "\n Increased patch number to $patchNumber \n"
git commit -m "Increased patch number to $patchNumber" python/GtBurst/version.py

echo "Generating file list with MD5s..."

#Now generate the list of files with their MD5, and upload it
find python -type f -name '*.py' -exec md5sum -b {} + >! __file_list
find data -type f ! -name '.*' -exec md5sum -b {} + >> __file_list

echo "Commiting the file list..."
git commit -m "Release of patch $patchNumber" __file_list

git push

cd ..

