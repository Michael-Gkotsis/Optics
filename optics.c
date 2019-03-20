#include <stdio.h>
#include <stdlib.h>
#include "FileHandling.h"
#include <math.h>
#include <string.h>
#include <time.h>


int main(int argc, char *argv[])
{

  int dim = 0;         // Dimensions of Elements
  int n = 0;           // n Elements
  int i,j,d;           // i counter for n, j counter for Neighborhood Elements, d counter for dimensions
  double minDist;      // Minimum Distance for a point to be assigned to clusters
  int minPoints;       // Minimum points for a cluster to be created
  int CounterNoise=0;  // Counter for Noise Elements
  int cluster = 0;     // Cluster index
  int h,k;             // h counter for Neighborhood Distance Calculation, k counter for files
  clock_t start, end;  // Variables for counting execution time
  char filename[40];
  char c;
  int choise = 0;
  double undefined = 9999;
  double min;
  int location;
  int size;
  int  e  = 0;
  double clustDist;      // Minimum Distance for a point to be assigned to clusters

/* -------------------------------------------------------------------------- */
    // Reading Dataset

   FILE* Dataset;

      printf("\n Give the DataSet file:");
      scanf("%s", filename);

    Dataset = fopen(filename, "r");
    if (!Dataset)
    {
        printf("\n There is something wrong with the Dataset file! \n\n");
        return -1;
    }

/* -------------------------------------------------------------------------- */





    dim = getColumns(Dataset); //Getting the dimensions of each Element of the Dataset
    rewind(Dataset);

    n = getRows(Dataset);       // Getting the Elements of the Dataset
    rewind(Dataset);

    printf("\n Elements:%d \n", n-1);
    printf("\n Dimensions:%d \n", dim);
     n--;
    printf("Give the amount of Minimum Points: ");
    scanf("%d",&minPoints );
    printf("Give the Generating Distance: ");
    scanf("%lf",&minDist);
    printf("Give the Clustering Distance: ");
    scanf("%lf",&clustDist);

/* -------------------------------------------------------------------------- */
        // All the necessary memory allocation

        float *X;   // Array of Elements
        X =(float *)calloc(n*dim, sizeof(float *));


        int *Cluster;
        Cluster =(int *)calloc(n,sizeof(int));


        int *visited; //Array for knowning which elements are visited and which are not
        visited =(int *) calloc(n,sizeof(int));



        float *distance; //array for holding distances of element i with each db element
        distance =(float *) calloc(n,sizeof(float));

        int *numPoints; //Array for holding the Neighborhood's size of each i element
        numPoints =(int *) calloc(n,sizeof(int));

        int *belong; //array for holding boolean value of an element belong to a Neighborhood
        belong =(int *) calloc(n,sizeof(int));

        int *nBelong; //array for holding boolean value of an element belong to a Neighborhood
        nBelong =(int *) calloc(n,sizeof(int));

        int *Noise; //Flag for points that are Noise
        Noise =(int *) calloc(n,sizeof(int));


        int *Core;  //Flag for core points,  0 for border, 1 for Core
        Core = (int *)calloc(n,sizeof(int));


        float *distance2 ;   //Array for holding Distances for each core point of a Neighborhood with Each Element of the dataset
        distance2 =(float *) calloc(n*n,sizeof(float));


        float *OrderList;   // Array of Elements
        OrderList =(float *) calloc(n*dim, sizeof(float));

        float *OrderReachabilityDistance; //Array for holding ordered reachability distances
        OrderReachabilityDistance = (float *)calloc(n,sizeof(float));

        float *coreDistance; //Array for holding core distances for each Core element
        coreDistance =(float *) calloc(n,sizeof(float));

        float *tempReachabilityDistance; //Array for holding ordered reachability distances
        tempReachabilityDistance =(float *) calloc(n,sizeof(float));

        float *reachabilityDistance; //Array for holding reachability distance for each element
        reachabilityDistance =(float *) calloc(n,sizeof(float));

        int *Seed;
        Seed =(int *) calloc(n,sizeof(int));

        for(i = n; i--;)
        {
        Cluster[i] = -1;
        reachabilityDistance[i] = undefined;
        coreDistance[i] = undefined;
        visited[i] = 0; //0 for unvisited , 1 for visited
        Seed[i] = 0;
        Noise[i] = 1;
        Core[i] = 0;
        }







/* -------------------------------------------------------------------------- */
                  // Passing elements to Array X[n][dim]



                    X = getData(Dataset,n,dim,X);



            fclose(Dataset);
  start = clock();
/* -------------------------------------------------------------------------- */
/* ---------------------------------OPTICS---------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
h = 0;
for(i = n; i--;)
{
  if(visited[i] != 1)
  {
   numPoints[i] = 0;
    visited[i] = 1;



for(d = dim; d--;)
  OrderList[h*dim + d] = X[i*dim + d];
  OrderReachabilityDistance[h] = undefined;



/* -------------------------------------------------------------------------- */
    //Finding the Neighborhood

       //Calculating the distance of Xi with each unvisited element
             for(j = n; j--;)
             {
               belong[j] = 0; // Mark temp 0 for current Point, used for reset
               Seed[j] = 0;
               distance[j] = 0; //Reset the distance before calculating for current point
               if(j != i)
               {


                    for(d = dim; d--;)
                    {
    //Calculating and storing distances of Xi element with the rest of the points
                   distance[j] += (X[j*dim + d] - X[i*dim + d])*(X[j*dim + d] - X[i*dim + d]);
                    }
                  distance[j] = sqrt(distance[j]);

              //Calculating the size of the Neighborhood
                      if(distance[j] <= minDist)
                      {
    //If Current elements belongs to the Neighborhood set temp 1 and Increase size
                         belong[j] = 1;
                            numPoints[i]++;
                      }
                }
              }
              // printf("%d : numPoints %d\n",i,numPoints[i] );
    //Raising by one so that the i element to be accounted for
       numPoints[i]++;

       if(numPoints[i] >= minPoints)
       {

        e = 0;

         for(j = n; j--;)
         {
           if(belong[j] != 0)
           {
     reachabilityDistance[j] = 0; //Reset the distance before calculating for current point

     for(d = dim; d--;)
     {
//Calculating and storing distances of Xi element with the rest of the points
    reachabilityDistance[j] += (X[j*dim + d] - X[i*dim + d])*(X[j*dim + d] - X[i*dim + d]);
     }
   reachabilityDistance[j] = sqrt(reachabilityDistance[j]);
  tempReachabilityDistance[e] = reachabilityDistance[j];
  e++;
  if(visited[j] != 1)
  Seed[j] = 1;
           }
         }

         qSort(tempReachabilityDistance,e);

       //  for(j = 0; j < e; j++)
       //  {
       //   printf("%d: ReachabilityDist : %lf \n",j,tempReachabilityDistance[j] );
       // }

       coreDistance[h] = tempReachabilityDistance[minPoints - 1];

       h++;

       // printf("%d CoreDist : %lf\n",i,coreDistance[i] );

       do{
          size = 0;

          min = 9999;
          for(j = n; j--;)
           {
               if(Seed[j] != 0)
               {
                 if(reachabilityDistance[j] < min)
                 {
                   min = reachabilityDistance[j];
                   location = j;
                 }
               }

           }



           Seed[location] = 0;
           visited[location] = 1;
           for(d = dim; d--;)
           OrderList[h*dim + d] = X[location*dim + d];

           OrderReachabilityDistance[h] = reachabilityDistance[location];
           numPoints[location] = 0;

           h++;

           for(k = n; k--;)
           {
             nBelong[k] = 0; // Mark temp 0 for current Point, used for reset
             distance[k] = 0; //Reset the distance before calculating for current point
             if(k != location)
             {

                  for(d = dim; d--;)
                  {
  //Calculating and storing distances of Xi element with the rest of the points
                 distance[k] += (X[k*dim + d] - X[location*dim + d])*(X[k*dim + d] - X[location*dim + d]);
                  }
                distance[k] = sqrt(distance[k]);

            //Calculating the size of the Neighborhood
                    if(distance[k] <= minDist)
                    {
  //If Current elements belongs to the Neighborhood set temp 1 and Increase size
                       nBelong[k] = 1;
                          numPoints[location]++;
                    }

                }
              }

            numPoints[location]++;

            if(numPoints[location] >= minPoints)
            {
              e = 0;

               for(k = n; k--;)
               {
                 if(nBelong[k] != 0)
                 {
                    if(Seed[k] != 0)
                    {
                      double temp = reachabilityDistance[k];

                      reachabilityDistance[k] = 0;

                      for(d = dim; d--;)
                      {
                 //Calculating and storing distances of Xi element with the rest of the points
                     reachabilityDistance[k] += (X[k*dim + d] - X[location*dim + d])*(X[k*dim + d] - X[location*dim + d]);
                      }
                    reachabilityDistance[k] = sqrt(reachabilityDistance[k]);

                    if(temp < reachabilityDistance[k])
                    {
                      reachabilityDistance[k] = temp;
                    }

                    tempReachabilityDistance[e] = reachabilityDistance[k];
                    e++;
                  }else
                  {
                    reachabilityDistance[k] = 0;

                    for(d = dim; d--;)
                    {
               //Calculating and storing distances of Xi element with the rest of the points
                   reachabilityDistance[k] += (X[k*dim + d] - X[location*dim + d])*(X[k*dim + d] - X[location*dim + d]);
                    }
                  reachabilityDistance[k] = sqrt(reachabilityDistance[k]);

                  tempReachabilityDistance[e] = reachabilityDistance[k];
                  e++;

                  if(visited[k] != 1)
                  Seed[k] = 1;
                  }
                 }
               }

               qSort(tempReachabilityDistance,e);
         coreDistance[h] = tempReachabilityDistance[minPoints - 1];

       }







          for(j = n; j--;)
          {
            if(Seed[j] == 1)
            {
              size++;
            }
          }
          // printf(" %d \n", size );
       }while(size != 0);



     }else
     {
       h++;
       coreDistance[h] = undefined;
     }
/* -------------------------------------------------------------------------- */




/* -------------------------------------------------------------------------- */
  }
}
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

cluster = 0;
for(j = h; j--; )
{
    if(OrderReachabilityDistance[j] > clustDist)
      {
        if(coreDistance[j] <= clustDist)
          {
            cluster++;
            Cluster[j] = cluster;
            Noise[j] = 0;

          }else
           {
             Noise[j] = 1;
           }
      }else
      {
            Cluster[j] = cluster;
            Noise[j] = 0;
      }
}


/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* ---------------------------------END-------------------------------------- */
/* -------------------------------------------------------------------------- */
end = clock();
/*----------------------------------------------------------------------------*/
double total_time = ((double) (end - start)) / CLOCKS_PER_SEC;



// for(i = 0; i < n ; i++)
// {
// printf("\n");
// printf("%d :",i+1 );
// printf("Noise: %d ",Noise[i]);
// printf("Core: %d ",Core[i]);
//  printf("Visited: %d ",visited[i]);
//  printf("Cluster %d ",Cluster[i] );
//    for(d = 0; d < dim; d++ )
//     printf("%lf ",X[i][d] );
//
//   }
CounterNoise = 0;
for(i = h; i--;)
{
  if(Noise[i] == 1)
  {
    CounterNoise++;
  }
}

printf("\n");
printf("Clusters : %d \n",cluster );
printf("\n");
printf("CounterNoise: %d\n",CounterNoise );
printf("\n");
printf("\n Time of Algorithm Execution: %lf \n\n",total_time);

printf("\n");

FILE* NoiseFile;

NoiseFile = fopen("FinalNoise.txt","w");

for(j = h; j--;){
  if(Noise[j] == 1)
  {
    for(d = dim; d--;)
    fprintf(NoiseFile, "%lf ",OrderList[j*dim + d] );

    fprintf(NoiseFile, "\n");
  }
}
fclose(NoiseFile);

FILE* ReachabilityDistFile;

ReachabilityDistFile = fopen("ReachabilityDist.txt","w");

for(j = h; j--;)
{
  fprintf(ReachabilityDistFile, "%lf \n",OrderReachabilityDistance[j]);
}

fclose(ReachabilityDistFile);

FILE* CoreFile;

CoreFile = fopen("CoreDist.txt","w");

for(j = h; j--;)
{
  fprintf(CoreFile, "%lf \n",coreDistance[j]);
}

fclose(CoreFile);

FILE* OrderFile;

OrderFile = fopen("OrderList.txt","w");

for(j = h; j--;)
{

  for(d = dim; d--;)
  fprintf(OrderFile, "%lf ",OrderList[j*dim + d]);

  fprintf(OrderFile, "\n");
}

fclose(OrderFile);


for(j = 200 ; j >= cluster ; j--)
{
char *fileName;
fileName = calloc(n,sizeof(char));
snprintf(fileName,n*sizeof(char),"C%d.txt",j);
remove(fileName);
free(fileName);
}

/* -------------------------------------------------------------------------- */

free(X);
free(Cluster);
free(visited);
free(Noise);
free(nBelong);
free(belong);
free(numPoints);
free(distance);
free(Core);
free(distance2);
free(OrderList);
free(OrderReachabilityDistance);
free(coreDistance);
free(tempReachabilityDistance);
free(reachabilityDistance);
free(Seed);


  return 0;
}
