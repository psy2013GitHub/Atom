
# Usage:
#       ~/sum_spec.sh -d l -o /media/dz/TOSHIBA/C/C_GFCD -i /media/dz/TOSHIBA/C/C_GFCD.nii
#!/bin/sh

# get arg
while getopts "i:o:d:" arg 
do
  case $arg in
       o)
        out_tlrc_img=$OPTARG
        ;;
       i)
        input_img=$OPTARG
        ;;
       d)
        direction=$OPTARG
        ;;
  esac
done


/usr/local/lib/AFNI/3dcopy $input_img $out_tlrc_img

/usr/local/lib/AFNI/3drefit -view tlrc -space MNI ${out_tlrc_img}+tlrc.HEAD

/usr/local/lib/AFNI/afni -niml &

cp ${out_tlrc_img}+tlrc.HEAD ${out_tlrc_img}+tlrc.BRIK ~/suma_mni/

echo $direction
if [ "$direction" = "l" ]; then
   /usr/local/lib/AFNI/suma -spec ~/suma_mni/N27_lh.spec -sv ~/suma_mni/MNI_N27+tlrc.HEAD
fi

if [ "$direction" = "r" ]; then
   /usr/local/lib/AFNI/suma -spec ~/suma_mni/N27_rh.spec -sv ~/suma_mni/MNI_N27+tlrc.HEAD
fi



