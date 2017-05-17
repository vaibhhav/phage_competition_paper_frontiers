/////////////////////////////////////////////////////////////////////////////
// Model of viral infection in bacterial colony
// Modified model of growth
// Bacteria don't move, grow at a fixed rate
// Each bacterium judges free available spots at random,
// and place progeny there (only if free spots are available)
// Bacteria collect damage every clock tick, and the mother retains damage
// upon binary fission. They die after after a threshold damage
// has been accumulated. Results are stored in "infectedResults.csv"
/////////////////////////////////////////////////////////////////////////////
#include <iostream>
// #include <conio.h>
#include <cmath>
// #include <windows.h>
#include <stdio.h>
#include <cstdlib>
#include <time.h>
#include <fstream>

using namespace std;

///////////////////////////////////////////////////////////////////////
// Parameters
//////////////////////////////////////////////////////////////////////
const int healthyGrowthRate = 20; //setting growth rate to one fission every 20 cycles of the internal clock of the bacterium
const int lysogenicGrowthRateA = 30; //7 //division time for infected bacteria
const int lysogenicGrowthRateB = 30;
const int healthyDeathTime = 150; //15 //healthy bacteria die every healthyDeathTime units
const int lysogenicDeathTimeA = 120; //22
const int lysogenicDeathTimeB = 120;
const int lysisBurstTimeA = 1000; //time (after infection) after which an infected lytic bacterium bursts
const int lysisBurstTimeB = 1000;
const float probInfectionA = 0.4; //6 //probability of infection per burst in bacterial vicinity
const float probInfectionB = 0.4;
const int decisionTime = 30; //time during which bacterium makes decision
const int birthTimeRange = 2; //+- values of time for which birth can occur
const int numberOfTriesA = 3;
const int numberOfTriesB = 3;
long bacteriaCount = 0; //counters keeping track of the numbers
const long latticeSize = 40;
long lysogenicBacteriaCountA = 0;
long lysogenicBacteriaCountB = 0;
long lyticBacteriaCountA = 0;
long lyticBacteriaCountB = 0;
long deathCount = 0; //counts number of dead bacteria
long maxTime = 200000;

///////////////////////////////////////////////////////////////////////
// Matrices
//////////////////////////////////////////////////////////////////////
long lattice[latticeSize][latticeSize] = {0}; //stores the status of each cell. 0: empty; 1: healthy; 2: lysogenic; 3: lytic
long infectionStatus[latticeSize][latticeSize] = {0}; //stores the infection counter for each healthy bacterium (see infectionThreshold)
long clockTicks[latticeSize][latticeSize]; //stores the age of each bacterium; acts as an internal clock
long lyticTimer[latticeSize][latticeSize];
long decisionState[latticeSize][latticeSize] = {0};
long multiplicityCounter[latticeSize][latticeSize]; //multiplicity of infection counter for each bacterium in decisionState 1
int allowed[4] = {0}; //keeps track of the free spots in the neighbourhood of each bacterium, where 0 is the top, and the numbers increase clockwise (in a plus shape)
int randomBox[4] = {-1, -1, -1, -1}; //array used to keep track of indices of allowed to which cell division is allowed
long healthyBirthTicker[latticeSize][latticeSize]; //keeps track of the random birth time around the division age for each healthy cell
long infectedBirthTicker[latticeSize][latticeSize]; //does the same as healthyBirthTicker, but for infected cells
ofstream results; //file storing the results
int p;

void reinitAllowed() //resets the value of allowed and randomBox every time a new bacterium is moved to
{
    for(int i = 0; i < 4; i++)
    {
        allowed[i] = 0;
        randomBox[i] = -1;
    }
}

///////////////////////////////////////////////////////////////////////
// Function prototypes
//////////////////////////////////////////////////////////////////////
void displayLattice(long);
void decideProgenySpot(int, long &, long &, long, long, int &);
int countFreeSpots(long &, long &);
//int countInfectors(long &, long &);
void generateBacterium(int, long, long);
void generateIsland(long, long, int);
void beginInfectionStageA(long, long);
void makeDecisionA(long, long);
void beginInfectionStageB(long, long);
void makeDecisionB(long, long);

int main()
{
    srand(10);
    //setting all clocks and tickers to null (-1)
    for(long x = 0; x < latticeSize; x++)
        for(long y = 0; y < latticeSize; y++)
        {
            clockTicks[x][y] = -1;
            healthyBirthTicker[x][y] = -1;
            infectedBirthTicker[x][y] = -1;
            lyticTimer[x][y] = -1;
        }

    //opening the results file and prepping it up
    results.open("infectedResults.csv", ios::out);
    results << "t,bacteriaCount,lysogenicBacteriaCountA,lyticBacteriaCountA,lysogenicBacteriaCountB,lyticBacteriaCountB,healthyBacteriaCount\n";

    ////////////////////////////////////////////////////////////////////////////////////
    // Initial conditions
    ////////////////////////////////////////////////////////////////////////////////////
    // 0 represents free spot, 1 represents healthy bacterium, 2 represents infected bacterium
    //generateIsland creates an island with an infected bacterium flanked by 4 healthy ones on an x shape centered at (x,y)
    generateIsland(latticeSize/4, latticeSize/4, 5);
    generateIsland(3*latticeSize/4, 3*latticeSize/4, 3);
    //generateIsland(latticeSize/4, 3*latticeSize/4, 1);
    //generateIsland(3*latticeSize/4, latticeSize/4, 1);

    for(long t = 0; t < maxTime || bacteriaCount == pow(latticeSize, 2); t++)
    {
        //srand(time(NULL));
        displayLattice(t); //Sleep(20000000000); //<-- I use sleep to pause the program and check some things randomly
        for(long x = 0; x < latticeSize; x++) //increment internal clocks for each bacterium, decrementing the birth tickers for each time step
            for(long y = 0; y < latticeSize; y++)
            {
                if(clockTicks[x][y] != -1)
                    clockTicks[x][y]++;
                if(healthyBirthTicker[x][y] >= 0)
                    healthyBirthTicker[x][y] -= 1;
                if(infectedBirthTicker[x][y] >= 0)
                    infectedBirthTicker[x][y] -= 1;
                if(lyticTimer[x][y] != -1)
                    lyticTimer[x][y]++;
            }

        //sweeping through the lattice
        for(long j = 0; j < latticeSize; j++)
        {for(long k = 0; k < latticeSize; k++)
            {
                if(lattice[j][k]) //looking for a non-empty site
                {
                    //making sure the ticker for the bacterium hasn't gone below -1
                    if(healthyBirthTicker[j][k] < -1)
                        healthyBirthTicker[j][k] = -1;
                    if(infectedBirthTicker[j][k] < -1)
                        infectedBirthTicker[j][k] = -1;

                    //this loop kills things if their time has come
                    if((clockTicks[j][k] == healthyDeathTime && lattice[j][k] == 1) || (clockTicks[j][k] == lysogenicDeathTimeA && lattice[j][k] == 2) || (clockTicks[j][k] == lysogenicDeathTimeB && lattice[j][k] == 4) || (lyticTimer[j][k] == lysisBurstTimeA && lattice[j][k] == 3) || (lyticTimer[j][k] == lysisBurstTimeB && lattice[j][k] == 5))
                    {
                        if(lattice[j][k] == 2)
                            lysogenicBacteriaCountA--;
                        else if(lattice[j][k] == 4)
                            lysogenicBacteriaCountB--;
                        else if(lattice[j][k] == 3)
                        {
                            lyticBacteriaCountA--;
                            beginInfectionStageA(j, k);
                        }
                        else if(lattice[j][k] == 5)
                        {
                            lyticBacteriaCountB--;
                            beginInfectionStageB(j, k);
                        }
                        lattice[j][k] = 0; //resetting all counters and stuff here
                        clockTicks[j][k] = -1;
                        healthyBirthTicker[j][k] = -1;
                        infectedBirthTicker[j][k] = -1;
                        multiplicityCounter[j][k] = 0;
                        infectionStatus[j][k] = 0;
                        decisionState[j][k] = 0;
                        bacteriaCount--;
                        deathCount++;
                        continue; //continue skips the rest of this iteration and begins the next one for this main for loop immediately
                    }

                    if(lattice[j][k] == 1)
                    {
                        if(infectionStatus[j][k] > 0) //negative infectionStatus implies infection from B phage, and positive, by A phage
                        {
                            p = rand() % 100;
                            if(p < (probInfectionA*100))
                            {
                                multiplicityCounter[j][k]++;
                                if(!decisionState[j][k])
                                    decisionState[j][k] = decisionTime;
                            }
                            infectionStatus[j][k]--;
                        }
                        else if(infectionStatus[j][k] < 0)
                        {
                            p = rand() % 100;
                            if(p < (probInfectionB*100))
                            {
                                multiplicityCounter[j][k]--;
                                if(!decisionState[j][k])
                                    decisionState[j][k] = -decisionTime;
                            }
                            infectionStatus[j][k]++;
                        }

                        if(decisionState[j][k] > 0) //again, positive decisionState implies decision via strategy A
                        {
                            decisionState[j][k]--;
                            if(!decisionState[j][k])
                                makeDecisionA(j, k);
                        }
                        else if(decisionState[j][k] < 0)
                        {
                            decisionState[j][k]++;
                            if(!decisionState[j][k])
                                makeDecisionB(j, k);
                        }
                    }

                    //choosing random times for cell division if the current clock of the cell reads a time just within the division age (so at 25 if the division age is 30 and the range is 5)
                    //srand(time(NULL)); //seeding the random number generator with the system time
                    if((healthyBirthTicker[j][k] == -1) && (clockTicks[j][k] != 0) && (clockTicks[j][k] % healthyGrowthRate == (healthyGrowthRate - birthTimeRange)) && lattice[j][k] == 1)
                        healthyBirthTicker[j][k] = rand() % (birthTimeRange*2 + 1); //rand() generates a random integer between 0 and whatever number is after the % sign (except that number itself)

                    //srand(time(NULL));
                    if((infectedBirthTicker[j][k] == -1) && (clockTicks[j][k] != 0) && (((clockTicks[j][k] % lysogenicGrowthRateA == (lysogenicGrowthRateA - birthTimeRange)) && lattice[j][k] == 2) || ((clockTicks[j][k] % lysogenicGrowthRateB == (lysogenicGrowthRateB - birthTimeRange)) && lattice[j][k] == 4)))
                        infectedBirthTicker[j][k] = rand() % (birthTimeRange*2 + 1);

                    //loop checks if the growth age for the cell has been reached and initiates cell division; clockTicks[j][k] != 0 makes sure newly-born bacteria can't pass this check
                    if(((healthyBirthTicker[j][k] == 0 && lattice[j][k] == 1) || (infectedBirthTicker[j][k] == 0 && (lattice[j][k] == 2 || lattice[j][k] == 4)) && clockTicks[j][k] != 0 && !infectionStatus[j][k] && !decisionState[j][k]))
                    {
                        long row, col; //keeps track of where the daughter should be placed
                        int flag = 0; //used in case of no free spots
                        int dice; //dice and tempDice are random variables that choose the division spot
                        int tempDice;
                        do
                        {
                            //srand(time(NULL)); //seeding the random number generator with the system time
                            if(countFreeSpots(j, k) == 0) //if no free spots
                            {   tempDice = -1; break;    } //break comes out of the loop
                            int boxCount = 0; //keeps the index of the traversal through randomBox
                            for(int index = 0; index < 4; index++)
                                if(allowed[index]) //allowed[index] is 1 if the index position is free (0 is top, 1 is right, 2 is bottom, 3 is left)
                                    randomBox[boxCount++] = index;
                            tempDice = randomBox[rand() % (boxCount)]; //picking out a random position from the free location indices stored in randomBox
                            break;
                        }while(tempDice != -1);
                        reinitAllowed();
                        dice = tempDice;

                        decideProgenySpot(dice, row, col, j, k, flag); //fixes values of row and col based on the position chosen by dice

                        if(!lattice[row][col] && !flag) //places progeny if free spots are available (!a is equivalent to a == 0)
                        {
                            lattice[row][col] = lattice[j][k];
                            clockTicks[row][col] = 0; //setting clock of the newborn to 0
                            bacteriaCount++;
                            healthyBirthTicker[row][col] = -1;
                            infectedBirthTicker[row][col] = -1;
                            if(lattice[row][col] == 2)
                                lysogenicBacteriaCountA++;
                            else if(lattice[row][col] == 4)
                                lysogenicBacteriaCountB++;
                            //displayLattice(t);
                        }
                    }
                }
            }
        }
    }
    cout << endl << "Simulation complete.";
    results.close();
    exit(0);
    // getch();
}

void generateBacterium(int type, long relativePosX, long relativePosY)
{
    lattice[latticeSize/2 + relativePosX][latticeSize/2 + relativePosY] = type;
    clockTicks[latticeSize/2 + relativePosX][latticeSize/2 + relativePosY] = 0;
    bacteriaCount++;
    if(type == 2)
        lysogenicBacteriaCountA++;
    if(type == 4)
        lysogenicBacteriaCountB++;
    if(type == 3)
    {
        lyticBacteriaCountA++;
        lyticTimer[latticeSize/2 + relativePosX][latticeSize/2 + relativePosY] = 0;
    }
    if(type == 5)
    {
        lyticBacteriaCountB++;
        lyticTimer[latticeSize/2 + relativePosX][latticeSize/2 + relativePosY] = 0;
    }
}

void generateIsland(long centerX, long centerY, int centerType)
{
    generateBacterium(centerType, centerX - latticeSize/2, centerY - latticeSize/2);
    generateBacterium(1, centerX - latticeSize/2 - 1, centerY - latticeSize/2 - 1);
    generateBacterium(1, centerX - latticeSize/2 - 1, centerY - latticeSize/2 + 1);
    generateBacterium(1, centerX - latticeSize/2 + 1, centerY - latticeSize/2 - 1);
    generateBacterium(1, centerX - latticeSize/2 + 1, centerY - latticeSize/2 + 1);
}

void displayLattice(long t) //prints the lattice to the screen and stores the data point in the results file
{
    /*system("cls");
    for(long y = 0; y < latticeSize; y++)
    {
        cout << endl;
        for(long x = 0; x < latticeSize; x++)
            cout << lattice[x][y] << " ";
    }*/

    results << t << "," << bacteriaCount << "," << lysogenicBacteriaCountA << "," << lyticBacteriaCountA << "," << lysogenicBacteriaCountB << "," << lyticBacteriaCountB << "," << (bacteriaCount-lysogenicBacteriaCountA-lyticBacteriaCountA-lysogenicBacteriaCountB-lyticBacteriaCountB) << "\n";
    /*cout << endl << endl << "Time: " << t;
    cout << endl << "Number of bacteria: " << bacteriaCount << " (" << (bacteriaCount - lysogenicBacteriaCount - lyticBacteriaCount) << " healthy + " << lysogenicBacteriaCount << " lysogenic + " << lyticBacteriaCount << " lytic.)";
    cout << endl << "Number of deaths in the colony: " << deathCount;*/
}

void decideProgenySpot(int dice, long &row, long &col, long j, long k, int &flag) //fixes row and col from dice
{
                        if(dice == -1)
                        {
                            flag = 1;
                        }
                        else if(dice == 0)
                        {
                            row = j-1;
                            col = k;
                        }
                        else if(dice == 1)
                        {
                            row = j;
                            col = k+1;
                        }
                        else if(dice == 2)
                        {
                            row = j+1;
                            col = k;
                        }
                        else if(dice == 3)
                        {
                            row = j;
                            col = k-1;
                        }
}

int countFreeSpots(long &j, long &k) //looks at all four neighbourhood spots and returns number of free spots and fills allowed with the status of each position
{
    int freeSpots = 0;
    if(!lattice[j-1][k] && j != 0)
    { freeSpots++; allowed[0] = 1; }
    if(!lattice[j+1][k] && j != latticeSize-1)
    { freeSpots++; allowed[2] = 1; }
    if(!lattice[j][k-1] && k != 0)
    { freeSpots++; allowed[3] = 1; }
    if(!lattice[j][k+1] && k != latticeSize-1)
    { freeSpots++; allowed[1] = 1; }
    return freeSpots;
}

void beginInfectionStageA(long j, long k)
{
    if(infectionStatus[j+1][k+1] == 0)
        infectionStatus[j+1][k+1] = numberOfTriesA;
    if(infectionStatus[j+1][k-1] == 0)
        infectionStatus[j+1][k-1] = numberOfTriesA;
    if(infectionStatus[j-1][k+1] == 0)
        infectionStatus[j-1][k+1] = numberOfTriesA;
    if(infectionStatus[j-1][k-1] == 0)
        infectionStatus[j-1][k-1] = numberOfTriesA;
    if(infectionStatus[j][k-1] == 0)
        infectionStatus[j][k-1] = numberOfTriesA;
    if(infectionStatus[j][k+1] == 0)
        infectionStatus[j][k+1] = numberOfTriesA;
    if(infectionStatus[j-1][k] == 0)
        infectionStatus[j-1][k] = numberOfTriesA;
    if(infectionStatus[j+1][k] == 0)
        infectionStatus[j+1][k] = numberOfTriesA;
}

void beginInfectionStageB(long j, long k)
{
    if(infectionStatus[j+1][k+1] == 0)
        infectionStatus[j+1][k+1] = -numberOfTriesB;
    if(infectionStatus[j+1][k-1] == 0)
        infectionStatus[j+1][k-1] = -numberOfTriesB;
    if(infectionStatus[j-1][k+1] == 0)
        infectionStatus[j-1][k+1] = -numberOfTriesB;
    if(infectionStatus[j-1][k-1] == 0)
        infectionStatus[j-1][k-1] = -numberOfTriesB;
    if(infectionStatus[j][k-1] == 0)
        infectionStatus[j][k-1] = -numberOfTriesB;
    if(infectionStatus[j][k+1] == 0)
        infectionStatus[j][k+1] = -numberOfTriesB;
    if(infectionStatus[j-1][k] == 0)
        infectionStatus[j-1][k] = -numberOfTriesB;
    if(infectionStatus[j+1][k] == 0)
        infectionStatus[j+1][k] = -numberOfTriesB;
}

void makeDecisionA(long j, long k)
{
    float probLysogeny = (exp(multiplicityCounter[j][k]) - 1)/(20*multiplicityCounter[j][k]);
    if(multiplicityCounter[j][k] > 3)
        probLysogeny = 1;
    if(multiplicityCounter[j][k] == 1)
        probLysogeny = 0.00;
    if(multiplicityCounter[j][k] == 2)
        probLysogeny = 1.00;
    if(multiplicityCounter[j][k] == 3)
        probLysogeny = 1.00;
    if(multiplicityCounter[j][k] > 3)
        probLysogeny = 1;
    if(rand() % 100 < (probLysogeny*100))
    {
        lattice[j][k] = 2;
        lysogenicBacteriaCountA++;
    }
    else
    {
        lattice[j][k] = 3;
        lyticTimer[j][k] = 0;
        lyticBacteriaCountA++;
    }
}

void makeDecisionB(long j, long k)
{
    multiplicityCounter[j][k] = abs(multiplicityCounter[j][k]);
    float probLysogeny = (exp(multiplicityCounter[j][k]) - 1)/(20*multiplicityCounter[j][k]);//0.4 * log(multiplicityCounter[j][k] + 1);
    if(multiplicityCounter[j][k] == 1)
        probLysogeny = 0.00;
    if(multiplicityCounter[j][k] == 2)
        probLysogeny = 0.00;
    if(multiplicityCounter[j][k] == 3)
        probLysogeny = 1.00;
    if(multiplicityCounter[j][k] > 3)
        probLysogeny = 1;
    if(rand() % 100 < (probLysogeny*100))
    {
        lattice[j][k] = 4;
        lysogenicBacteriaCountB++;
    }
    else
    {
        lattice[j][k] = 5;
        lyticTimer[j][k] = 0;
        lyticBacteriaCountB++;
    }
}
