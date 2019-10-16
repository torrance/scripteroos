#! /bin/bash
set -x
set -e

dir=$(dirname $1)

rm -r $dir/tmp.im $dir/tmp-regridded.im $dir/tmp.fits $dir/tmp-template.fits || true
ln -s $(basename $1) $dir/tmp.fits

ext=${2##*.}
if [[ $ext != 'im' ]]; then
  echo "Coverting template file..."
  ln -s $(realpath $2) $dir/tmp-template.fits
  $MIRDIR/fits in=${dir}/tmp-template.fits out=${dir}/tmp-template.im op=xyin
else
  echo "Linking to template file..."
  ln -s $(realpath $2) $dir/tmp-template.im
fi

$MIRDIR/fits in=${dir}/tmp.fits out=${dir}/tmp.im op=xyin
$MIRDIR/regrid in=${dir}/tmp.im out=${dir}/tmp-regridded.im tin=${dir}/tmp-template.im
$MIRDIR/fits in=${dir}/tmp-regridded.im out=${dir}/tmp-regridded.fits op=xyout
mv ${dir}/tmp-regridded.fits $(dirname $1)/$(basename $1 .fits)-regridded.fits

rm -r ${dir}/tmp.im ${dir}/tmp-regridded.im ${dir}/tmp.fits ${dir}/tmp-template.fits ${dir}/tmp-template.im || true

