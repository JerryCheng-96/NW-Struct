using namespace std;

extern "C" {
    int Blosum_NW_Relay(char*, int, char*, int, float(*sim_func)(void*, void*), int, int);
}

int Blosum_NW_Relay(char* chr_seq_1, int len_seq_1, char* chr_seq_2, int len_seq_2,
                    float (*sim_func)(void*, void*), int gap_open, int gap_ext) {
    void** seq_1 = (void**)malloc(sizeof(void*) * (len_seq_1 + 10));
    void** seq_2 = (void**)malloc(sizeof(void*) * (len_seq_2 + 10));

    for (int i = 0; i < len_seq_1; i++) {
        seq_1[i] = (void*)(chr_seq_1 + i);
    }

    for (int i = 0; i < len_seq_2; i++) {
        seq_2[i] = (void*)(chr_seq_2 + i);
    }

    NW_Align(seq_1, len_seq_1, seq_2, len_seq_2, sim_func, gap_open, gap_ext);
}