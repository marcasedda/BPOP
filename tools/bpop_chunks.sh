i=0
N=4
while [ "$i" -lt "$N" ];
do
  echo "$i" > chunk.out
  echo "cp BPOP.x chnk"$i""
  mkdir chnk"$i"
  cp BPOP.x chunk.out chnk"$i"
  cd chnk"$i"
  ./BPOP.x > Chunk_"$i"_Output.txt &
  cd ../
  i=$(($i+1))

done
