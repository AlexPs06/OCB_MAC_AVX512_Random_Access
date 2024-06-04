#include <iostream>
#include <wmmintrin.h>
#include <immintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>

#define ALIGN(n) __attribute__ ((aligned(n)))
#define pipeline 1

#define EXPAND_ASSIST(v1,v2,v3,v4,shuff_const,aes_const)                    \
    v2 = _mm_aeskeygenassist_si128(v4,aes_const);                           \
    v3 = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(v3),              \
                                         _mm_castsi128_ps(v1), 16));        \
    v1 = _mm_xor_si128(v1,v3);                                              \
    v3 = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(v3),              \
                                         _mm_castsi128_ps(v1), 140));       \
    v1 = _mm_xor_si128(v1,v3);                                              \
    v2 = _mm_shuffle_epi32(v2,shuff_const);                                 \
    v1 = _mm_xor_si128(v1,v2)

using namespace std;

static void PMAC(unsigned char *K_1, unsigned char *N,unsigned char *M, int size, unsigned char *T);
static void AES_128_Key_Expansion(const unsigned char *userkey, void *key);
static inline void AES_encrypt(__m128i tmp, __m128i *out,__m128i *key, unsigned rounds);
static void AES_cast_128_to_512_key2(__m128i *key,__m512i *key_512);
static void imprimiArreglo(int tam, unsigned char *in );
static inline void AES_encrypt_512(__m512i tmp, __m512i *out,__m512i *key, unsigned rounds);

// int main(){

//     ALIGN(64) unsigned char plaintext[128]=  {
//                                              0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
//                                              0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
//                                              0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
//                                              0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
                                             
//                                              0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
//                                              0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
//                                              0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
//                                              0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 

//                                             };
//     ALIGN(16) unsigned char tag[16 ]={ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
//     ALIGN(16) unsigned char K_1[16 ]={ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
//     ALIGN(16) unsigned char N[16 ]={ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};

//     PMAC(K_1, N, plaintext, 64, tag);
    
//     printf("\n");
//     imprimiArreglo(16, tag);
//     return 0;
// }


 char infoString[]= "PMAC AVX512 Random Access";  /* Each AE implementation must have a global one */

#ifndef MAX_ITER
#define MAX_ITER 16384
#endif

int main(int argc, char **argv)
{
	/* Allocate locals */
	ALIGN(64) unsigned char pt[MAX_ITER] = {0};
	ALIGN(16) unsigned char key[16]={ 0x00,0x01,0x02,0x03,
                                       0x04,0x05,0x06,0x07,
                                       0x08,0x09,0x0a,0x0b,
                                       0x0c,0x0d,0x0e,0x0f};

     ALIGN(16) unsigned char tag[16 ]={ 0x00,0x01,0x02,0x03,
                                       0x04,0x05,0x06,0x07,
                                       0x08,0x09,0x0a,0x0b,
                                       0x0c,0x0d,0x0e,0x0f};
    ALIGN(16) unsigned char nonce[16 ]={ 0x00,0x01,0x02,0x03,
                                       0x04,0x05,0x06,0x07,
                                       0x08,0x09,0x0a,0x0b,
                                       0x0c,0x0d,0x0e,0x0f};
	char outbuf[MAX_ITER*15+1024];
	int iter_list[2048]; /* Populate w/ test lengths, -1 terminated */
	char *outp = outbuf;
	int iters, i, j, len;
	double Hz,sec;
	double ipi=0, tmpd;
	clock_t c;
	iter_list[0] = 64;
	iter_list[1] = 128;
	iter_list[2] = 256;
	iter_list[3] = 512;
	iter_list[4] = 1024;
	iter_list[5] = 2048;
	iter_list[6] = 4096;
	iter_list[7] = 8192;
	iter_list[8] = 16384;
	iter_list[9] = -1;
    /* Create file for writing data */
	FILE *fp = NULL;
    char str_time[25];
	time_t tmp_time = time(NULL);
	struct tm *tp = localtime(&tmp_time);
	strftime(str_time, sizeof(str_time), "%F %R", tp);
	if ((argc < 2) || (argc > 3)) {
		printf("Usage: %s MHz [output_filename]\n", argv[0]);
		return 0;
	} else {
		Hz = 1e6 * strtol(argv[1], (char **)NULL, 10);
		if (argc == 3)
			fp = fopen(argv[2], "w");
	}

	
    outp += sprintf(outp, "%s ", infoString);

    outp += sprintf(outp, "- Run %s\n\n",str_time);

	// outp += sprintf(outp, "Context: %d bytes\n");
    
	printf("Starting run...\n");fflush(stdout);


	// /*
	//  * Get time for key setup
	//  */
	// iters = (int)(Hz/520);
	// do {
	
	// 	c = clock();
	// 	for (j = 0; j < iters; j++) {

	// 	}
	// 	c = clock() - c;
	// 	sec = c/(double)CLOCKS_PER_SEC;
	// 	tmpd = (sec * Hz) / (iters);
		
	// 	if ((sec < 1.2)||(sec > 1.3))
	// 		iters = (int)(iters * 5.0/(4.0 * sec));
	// 	printf("%f\n", sec);
	// } while ((sec < 1.2) || (sec > 1.3));

	
	printf("key -- %.2f (%d cycles)\n",sec,(int)tmpd);fflush(stdout);
	outp += sprintf(outp, "Key setup: %d cycles\n\n", (int)tmpd);

	/*
	 * Get times over different lengths
	 */
	iters = (int)(Hz/1000);
	i=0;
	len = iter_list[0];
	while (len >= 0) {
	
		do {
		

            PMAC(key,nonce,pt,iter_list[i],tag);

			c = clock();
			for (j = 0; j < iters; j++) {
                PMAC(key,nonce,pt,iter_list[i],tag);
			}
			c = clock() - c;
			sec = c/(double)CLOCKS_PER_SEC;
			tmpd = (sec * Hz) / ((double)len * iters);
			
			if ((sec < 1.2)||(sec > 1.3))
				iters = (int)(iters * 5.0/(4.0 * sec));
			
		} while ((sec < 1.2) || (sec > 1.3));
		
		printf("%d -- %.2f  (%6.2f cpb)\n",len,sec,tmpd);fflush(stdout);
		outp += sprintf(outp, "%5d  %6.2f\n", len, tmpd);
		if (len==44) {
			ipi += 0.05 * tmpd;
		} else if (len==552) {
			ipi += 0.15 * tmpd;
		} else if (len==576) {
			ipi += 0.2 * tmpd;
		} else if (len==1500) {
			ipi += 0.6 * tmpd;
		}
		
		++i;
		len = iter_list[i];
	}
	outp += sprintf(outp, "ipi %.2f\n", ipi);
	if (fp) {
        fprintf(fp, "%s", outbuf);
        fclose(fp);
    } else
        fprintf(stdout, "%s", outbuf);

	return ((pt[0]==12) && (pt[10]==34) && (pt[20]==56) && (pt[30]==78));
}

static void PMAC(unsigned char *K_1, unsigned char *N,unsigned char *M, int size, unsigned char *T){

    int m_blocks = 0;
    if (size%64==0)
        m_blocks=(size/64);
    else
        m_blocks=(size/64) + 1;

    static __m512i * plain_text = (__m512i*) M;
   
    __m128i nonce = _mm_setzero_si128();
    __m512i nonce_512;
    __m512i nonce_temp[1];
    __m128i Tag = _mm_setzero_si128();;
    __m128i keys_128[11];
    __m512i keys_512[11];
    __m512i keys_0[3];
    __m128i keys_0_128[3];

    __m512i S_temp;

    __m512i sum_nonce = _mm512_set_epi64(0,4, 0,4, 0,4, 0,4);

    union {__m128i bl128[4]; __m512i bl512;} S;

    keys_0[0] = _mm512_setzero_si512();
    keys_0[1] = _mm512_setzero_si512();
    keys_0[2] = _mm512_setzero_si512();
    keys_0_128[0] = _mm_setzero_si128();
    keys_0_128[1] = _mm_setzero_si128();
    keys_0_128[2] = _mm_setzero_si128();


    AES_128_Key_Expansion(K_1, keys_128);
    AES_cast_128_to_512_key2(keys_128, keys_512);

    nonce = _mm_setzero_si128();
    nonce = _mm_load_si128((__m128i *)&N[0]);
    
    AES_encrypt(nonce, &nonce, keys_128, 10);

    for (size_t i = 0; i < 4; i++){
        S.bl128[i]=nonce;
    }
    
    sum_nonce = _mm512_set_epi64(0,3, 0,2, 0,1, 0,0);

    nonce_512 = _mm512_add_epi64(sum_nonce, S.bl512);

    S.bl512 = _mm512_setzero_si512();
    sum_nonce = _mm512_set_epi64(0,4, 0,4, 0,4, 0,4);

    for (size_t i = 0; i < m_blocks; i++){

        nonce_temp[0]=nonce_512; 
        
        AES_encrypt_512(nonce_temp[0], &nonce_temp[0], keys_0, 2);
        
        plain_text[i] =_mm512_xor_si512(plain_text[i],nonce_temp[0]);
        
        AES_encrypt_512(plain_text[i], &plain_text[i], keys_512, 10);

        S_temp=_mm512_xor_epi64(plain_text[i],S_temp);        
        nonce_512=_mm512_add_epi64(nonce_512, sum_nonce);

    }


    S.bl512=S_temp;
    for (size_t i = 0; i < 4; i++){
        Tag=_mm_xor_si128(Tag,S.bl128[i]);
    }
    S.bl512=nonce_512;
    nonce = S.bl128[0];

    AES_encrypt(nonce, &nonce, keys_0_128, 2);

    AES_encrypt(nonce, &nonce, keys_128, 10);
    
    Tag=_mm_xor_si128(Tag,nonce);
	_mm_store_si128 ((__m128i*)T,Tag);
}


static inline void AES_encrypt(__m128i tmp, __m128i *out,__m128i *key, unsigned rounds){
	int j;
	tmp = _mm_xor_si128 (tmp,key[0]);
	for (j=1; j<rounds; j++)  tmp = _mm_aesenc_si128 (tmp,key[j]);
	tmp = _mm_aesenclast_si128 (tmp,key[j]);
	_mm_store_si128 ((__m128i*)out,tmp);
}


static inline void AES_encrypt_512(__m512i tmp, __m512i *out,__m512i *key, unsigned rounds){
	int j;
	tmp = _mm512_xor_si512 (tmp,key[0]);
	for (j=1; j<rounds; j++)  tmp = _mm512_aesenc_epi128 (tmp,key[j]);
	tmp = _mm512_aesenclast_epi128 (tmp,key[j]);
	_mm512_store_si512 ((__m128i*)out,tmp);
}



static void AES_128_Key_Expansion(const unsigned char *userkey, void *key)
{
    __m128i x0,x1,x2;
    __m128i *kp = (__m128i *)key;
    kp[0] = x0 = _mm_loadu_si128((__m128i*)userkey);
    x2 = _mm_setzero_si128();
    EXPAND_ASSIST(x0,x1,x2,x0,255,1);   kp[1]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,2);   kp[2]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,4);   kp[3]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,8);   kp[4]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,16);  kp[5]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,32);  kp[6]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,64);  kp[7]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,128); kp[8]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,27);  kp[9]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,54);  kp[10] = x0;
}
static void AES_cast_128_to_512_key2(__m128i *key,__m512i *key_512){
    union {__m128i oa128[4]; __m512i oa512;} oa;
    for(int i = 0; i< 11; i++ ){
        oa.oa128[0] = key[i];
        oa.oa128[1] = key[i];
        oa.oa128[2] = key[i];
        oa.oa128[3] = key[i];
        key_512[i]=oa.oa512;
    }

}


void imprimiArreglo(int tam, unsigned char *in )
{

    for (int i = 0; i<tam; i++){
        printf("%02x", in[i] );
    }
    printf("\n" );

}