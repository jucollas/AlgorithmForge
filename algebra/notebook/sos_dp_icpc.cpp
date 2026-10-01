//SUBSET convolution
for(int e=1;e<n;e*=2)for(int i=0;i<n;++i)
  if(i&e)vec[i  ]+=tp?vec[i^e]:-vec[i^e];
//SUPERSET convolution
for(int e=1;e<n;e*=2)for(int i=0;i<n;++i)
  if(i&e)vec[i^e]+=tp?vec[i  ]:-vec[i  ];
