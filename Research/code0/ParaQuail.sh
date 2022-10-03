#!/bin/sh

#HELLO TEST!

rm -r ../extract/Mesh$JIND

#Build x points dependent on node
xpoints[$JIND]=$((($PROG*8+$JIND)*4))

echo ${xpoints[*]}

#Build ypoints, same for all nodes


y=0
Yf=63
while [ $y -le $Yf ] 
do
  ypoints[$y]=$(($y*4))
  y=$(($y+1))
done
echo ${ypoints[*]}

#Manipulate periodic.f
x=0
n=0
char1=9  # Splits Cases0-9 and Cases10-32  


mkdir ../extract/Mesh$JIND
sed -i -e '1179s/.*/      'si=${xpoints[JIND]}'/' periodic.f
while [ $n -le $Yf ]
do
  sed -i -e '1181s/.*/      'sk=${ypoints[n]}'/' periodic.f
  make
  cd ../extract/Mesh$JIND/
  mkdir Case$n
  cd Case$n
  cp ../../../code$JIND/case1/input.dat .
  cp ../../../code$JIND/case1/input_per.dat .
  cp ../../../code$JIND/case1/diablo.start .
  ../../../code$JIND/diablo

    cd ../../
    if [ $n -le $char1 ]
        then
        sed -i "24s/.*/      character*11 cases/" res_to_dat.f
    else
        sed -i "24s/.*/      character*12 cases/" res_to_dat.f
    fi
    sed -i "57s/.*/      cases='Mesh$JIND\/Case$n'/" res_to_dat.f
    make
    ./res_to_dat
    cd Mesh$JIND/Case$n
    fil="Th01250"
    if [ -e "$fil" ]
	then
        rm dia*
    fi
  cd ../../../code$JIND/
  n=$(($n+1))
done


