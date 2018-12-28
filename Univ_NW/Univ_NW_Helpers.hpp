using namespace std;

void  print_matrix( float ** F, char* seq_1, char* seq_2, int len_seq_1, int len_seq_2 );
void  print_int_matrix( int ** F, char* seq_1, char* seq_2, int len_seq_1, int len_seq_2 );

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


void  print_int_matrix( int ** F, char* seq_1, char* seq_2, int len_seq_1, int len_seq_2 )
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

