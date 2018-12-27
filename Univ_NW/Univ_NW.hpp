using namespace std;

extern "C" {
    void  print_matrix( float ** F, char* seq_1, char* seq_2, int len_seq_1, int len_seq_2 );
    int arg_max(float* theArray, float* pMaxValue, int len);
    void NW_Align(void** seq_1, char* chr_seq_1, int len_seq_1, void** seq_2, char* chr_seq_2, int len_seq_2, float (*sim_func)(void*, void*), int gap_start, int gap_ext);
    void ReverseString(char * theString);
}


void  print_matrix( float ** F, char* seq_1, char* seq_2, int len_seq_1, int len_seq_2 )
{
        int  L1 = len_seq_1;
        int  L2 = len_seq_2;

        cout << "        ";
        for( int j = 0; j < L1; j++ )
        {
                cout << seq_1[ j ] << "   ";
        }
        cout << "\n  ";

        for( int i = 0; i <= L2; i++ )
        {
                if( i > 0 )
                {
                        cout << seq_2[ i-1 ] << " ";
                }
                for( int j = 0; j <= L1; j++ )
                {
                        cout.width( 3 );
                        cout << F[ i ][ j ] << " ";
                }
                cout << endl;
        }
}



int arg_max(float* theArray, float* pMaxValue, int len) {
    float maxValue = theArray[0];
    int maxIndex = 0;

    for(int i = 0; i < len; i++) {
        if (theArray[i] > maxValue) {
            maxValue = theArray[i];
            maxIndex = i;
        }
    }

    *pMaxValue = maxValue;
    return maxIndex;
}

void ReverseString(char * theString)
{
	char temp = '\0';
	int length = strlen(theString);

	for (int i = 0; i < length / 2 - 1; i++)
	{
		temp = theString[i];
		theString[i] = theString[length - 1 - i];
		theString[length - 1 - i] = temp;
	}
}

void NW_Align(void** seq_1, char* chr_seq_1, int len_seq_1, 
              void** seq_2, char* chr_seq_2, int len_seq_2, 
              float (*sim_func)(void*, void*), 
              int gap_start, int gap_ext) {
    // Printing parameters
    printf("%s\n%s\n%d\t%d\n%d\t%d\n", chr_seq_1, chr_seq_2, len_seq_1, len_seq_2, gap_start, gap_ext);


    // Initializing the matrices
    // The rows / first layer pointers
    float** scoresMat = (float**)malloc((len_seq_2 + 2) * sizeof(float*));
    int** dirMat = (int**)malloc((len_seq_2 + 2) * sizeof(int*));

    // The cols
    scoresMat[0] = (float*)malloc((len_seq_1 + 2) * (len_seq_2 + 2) * sizeof(float));
    dirMat[0] = (int*)malloc((len_seq_1 + 2) * (len_seq_2 + 2) * sizeof(int));
    printf("Size of the matrices: %d\n", (len_seq_1 + 2) * (len_seq_2 + 2));
    for (int row = 1; row < len_seq_2 + 2; row++) {
        scoresMat[row] = scoresMat[0] + row * (len_seq_1 + 2);
        dirMat[row] = dirMat[0] + row * (len_seq_1 + 2);
    }

    //for(int i = 0; i <= len_seq_1 + 1; i++) {
    //    scoresMat[i] = (float*)malloc((len_seq_2 + 2) * sizeof(float));
    //    dirMat[i] = (char*)malloc((len_seq_2 + 2) * sizeof(char));
    //}


    // Setting init values
    scoresMat[0][0] = 0;
    scoresMat[0][1] = -gap_start;
    scoresMat[1][0] = -gap_start;
    dirMat[0][0] = 0;
    dirMat[0][1] = 1;
    dirMat[1][0] = -1;

    for(int j = 2; j <= len_seq_1; j++) {
        scoresMat[0][j] = -(gap_ext * (j - 1)) - gap_start;
        dirMat[0][j] = j;
    }
    
    for(int i = 2; i <= len_seq_2; i++) {
        scoresMat[i][0] = -(gap_ext * (i - 1)) - gap_start;
        dirMat[i][0] = -i;
    }

    float scoreValues[] = {0,0,0};

    for(int i = 1; i <= len_seq_2; i++) {
        for(int j = 1; j <= len_seq_1; j++) {
            // Considering "local" (2x2) elements
            float sim_val = sim_func(seq_2[i - 1], seq_1[j - 1]);
            scoreValues[0] = scoresMat[i][j - 1] - gap_start;
            scoreValues[1] = scoresMat[i - 1][j - 1] + sim_val;
            scoreValues[2] = scoresMat[i - 1][j] - gap_start;

            float max_argval = 0.0;
            int max_arg = arg_max(scoreValues, &max_argval, 3);
            scoresMat[i][j] = max_argval;
            dirMat[i][j] = max_arg - 1;

            // Considering "distance" (row/col) elements
            // A distant, previous **ROW** elem "travels" here?
            for (int k = 0; k < j; k++) {
                float theScore = scoresMat[i][k] - (j-k-1) * gap_ext - gap_start;
                if (theScore > scoresMat[i][j]) {
                    scoresMat[i][j] = theScore;
                    dirMat[i][j] = j - k;
                }
            }

            // A distant, previous **COL** elem "travels" here?
            for (int k = 0; k < i; k++) {
                float theScore = scoresMat[k][j] - (i-k-1) * gap_ext - gap_start;
                if (theScore > scoresMat[i][j]) {
                    scoresMat[i][j] = theScore;
                    dirMat[i][j] = k - i;
                }
            }

            //printf("Now at (%d, %d)\n", i, j);
            //printf("Now at (%d, %d), scoreValues = (%f, %f, %f), max_arg = %d, max_val = %f. \n", i, j, scoreValues[0], scoreValues[1], scoreValues[2], max_arg, max_argval);
            //printf("\tD: %f, %d; L: %f, %d; T: %f, %d. SimValue(%c, %c) = %f\n", scoresMat[i-1][j-1], dirMat[i-1][j-1], scoresMat[i][j-1], dirMat[i][j-1], scoresMat[i-1][j], dirMat[i-1][j], chr_seq_2[i-1], chr_seq_1[j-1], sim_val);
        }
    }



    print_matrix(scoresMat, chr_seq_1, chr_seq_2, len_seq_1, len_seq_2);
    return;




    // Finding start point

    int max_row = 0;
    int max_col = len_seq_2;
    int max_val = scoresMat[0][len_seq_2];
    
    for(int i = 1; i <= len_seq_1; i++) {
        if (scoresMat[i][len_seq_2] > max_val) {
            max_row = i;
            max_col = len_seq_2;
            max_val = scoresMat[i][len_seq_2];
        }
    }

    for(int j = 1; j <= len_seq_2; j++) {
        if (scoresMat[len_seq_1][j] > max_val) {
            max_row = len_seq_1;
            max_col = j;
            max_val = scoresMat[len_seq_1][j];
        }
    }

    
    // Finding the path
    char* seq_1_res = (char*)malloc(sizeof(char) * (len_seq_1 + len_seq_2));
    char* seq_2_res = (char*)malloc(sizeof(char) * (len_seq_1 + len_seq_2));
    int idx_row = max_row;
    int idx_col = max_col;
    int idx_1_res = 0;
    int idx_2_res = 0;
    int idx_1_seq = len_seq_1 - 1;
    int idx_2_seq = len_seq_2 - 1;

    for (int i  = 0; i < len_seq_1 + len_seq_2; i++) {
        seq_1_res[i] = '\0';
        seq_2_res[i] = '\0';
    }

    if (max_col == len_seq_2) {     // The 2nd seq 'begins backwards' (ends) with gaps
        for (int i = 0; i < len_seq_1 - max_row; i++) {
            seq_1_res[i] = chr_seq_1[idx_1_seq];
            seq_2_res[i] = '-';
            idx_1_seq--;
        }
        idx_1_res = len_seq_1 - max_row;
        idx_2_res = len_seq_1 - max_row;
    }

    if (max_row == len_seq_1) {     // The 1st seq 'begins backwards' (ends) with gaps
        for (int i = 0; i < len_seq_2 - max_col; i++) {
            seq_1_res[i] = '-';
            seq_2_res[i] = chr_seq_2[idx_2_seq];
            idx_2_seq--;
        }
        idx_1_res = len_seq_2 - max_col;
        idx_2_res = len_seq_2 - max_col;
    }

    while (idx_row != 0 || idx_col != 0) {
        switch (dirMat[idx_row][idx_col])
        {
            case 'd':   // Both steps
                seq_1_res[idx_1_res] = chr_seq_1[idx_1_seq];
                seq_2_res[idx_2_res] = chr_seq_2[idx_2_seq];
                idx_1_res++;
                idx_2_res++;
                idx_1_seq--;
                idx_2_seq--;
                idx_col--;
                idx_row--;
                break;
            
            case 'l':   // Only seq_2 steps
                seq_1_res[idx_1_res] = '-';
                seq_2_res[idx_2_res] = chr_seq_2[idx_2_seq];
                idx_1_res++;
                idx_2_res++;
                idx_2_seq--;
                idx_col--;
                break;

             case 't':  // Only seq_1 steps
                seq_1_res[idx_1_res] = chr_seq_1[idx_1_seq];
                seq_2_res[idx_2_res] = '-';
                idx_1_res++;
                idx_2_res++;
                idx_1_seq--;
                idx_row--;
                break;
        
            default:
                break;
        }
    }

    ReverseString(seq_1_res);
    ReverseString(seq_2_res);

    // Output

    printf("\n\n%s\n", seq_1_res);
    printf("%s\n", seq_2_res);

    return;
}