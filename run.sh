#RED='\033[0;31m';
#NC='\033[0m';
#NL='\n';

if sudo g++ -fopenmp main.cpp algos.cpp -o program;
then
  echo "compilation was successful"

  if sudo ./program;
  then
    echo "run successful"
  else
    echo "runtime error"
  fi

else
  echo "compilation problems"
fi

sudo chmod a+rwx ./*
