
dirout=/mnt/netapp-barrier/public_html/sphinx/source/

cp -prfv _static/*png _static/*py ${dirout}_static/
cp -prfv conf.py ${dirout}

for f in *rst
do
    fout=`basename $f`
    if [ "$fout" != "index.rst" ]
    then
        cp -prvf $f ${dirout} 
    fi

done

