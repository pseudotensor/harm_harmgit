#!/bin/sh
#prefa=`echo $fil | sed "s/\./ /"`
#pref=`echo $prefa | awk '{print $1}'`
oldname=$1
newname=`echo $oldname | sed "s/def/dec/"`
echo $oldname to $newname
sed 's/^PFTYPE /extern PFTYPE /' $oldname > decs.${oldname}.temp11.h
sed 's/^SFTYPE /extern SFTYPE /' decs.${oldname}.temp11.h > decs.${oldname}.temp1.h
sed 's/^FTYPE /extern FTYPE /' decs.${oldname}.temp1.h > decs.${oldname}.temp1half.h
sed 's/^void /extern void /' decs.${oldname}.temp1half.h > decs.${oldname}.temp1half2.h
sed 's/^CTYPE /extern CTYPE /' decs.${oldname}.temp1half2.h > decs.${oldname}.temp2half.h
sed 's/^struct /extern struct /' decs.${oldname}.temp2half.h > decs.${oldname}.temp2.h
sed 's/^extern Sextern FTYPE /extern SFTYPE /' decs.${oldname}.temp2.h > decs.${oldname}.temp22.h
sed 's/^extern Pextern FTYPE /extern PFTYPE /' decs.${oldname}.temp22.h > decs.${oldname}.temp3.h
sed 's/^int /extern int /' decs.${oldname}.temp3.h > decs.${oldname}.temp4.h
sed 's/^short /extern short /' decs.${oldname}.temp4.h > decs.${oldname}.temp5.h
sed 's/^MPI_Status /extern MPI_Status /' decs.${oldname}.temp5.h > decs.${oldname}.temp6.h
sed 's/^MPI_Group /extern MPI_Group /' decs.${oldname}.temp6.h > decs.${oldname}.temp7.h
sed 's/^MPI_Comm /extern MPI_Comm /' decs.${oldname}.temp7.h > decs.${oldname}.temp8.h
sed 's/^char /extern char /' decs.${oldname}.temp8.h > decs.${oldname}.temp9.h
sed 's/^FILE\*/extern FILE\* /' decs.${oldname}.temp9.h > decs.${oldname}.temp92.h
sed 's/FILE \*/extern FILE\* /' decs.${oldname}.temp92.h > decs.${oldname}.temp10.h
sed 's/^mpidefs.h/mpidecs.h/' decs.${oldname}.temp10.h > decs.${oldname}.temp11.h
sed 's/^double /extern double /' decs.${oldname}.temp11.h > decs.${oldname}.temp12.h
sed 's/^float /extern float /' decs.${oldname}.temp12.h > decs.${oldname}.temp13.h
sed 's/^long /extern long /' decs.${oldname}.temp13.h > decs.${oldname}.temp14.h
sed 's/defs/decs/' decs.${oldname}.temp14.h > $newname
rm decs.${oldname}.temp*.h
