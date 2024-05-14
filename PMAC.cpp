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

int main(){

    ALIGN(64) unsigned char plaintext[128]=  {
                                             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 
                                             
                                             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
                                             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 

                                            };
    ALIGN(16) unsigned char tag[16 ]={ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
    ALIGN(16) unsigned char K_1[16 ]={ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
    ALIGN(16) unsigned char N[16 ]={ 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};

    PMAC(K_1, N, plaintext, 64, tag);
    
    printf("\n");
    imprimiArreglo(16, tag);
    return 0;
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
    __m512i keys_0[2];
    __m512i S_temp;

    __m512i sum_nonce = _mm512_set_epi64(0,4, 0,4, 0,4, 0,4);

    union {__m128i bl128[4]; __m512i bl512;} S;

    keys_0[0] = _mm512_setzero_si512();
    keys_0[1] = _mm512_setzero_si512();


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



    // imprimiArreglo(16,(unsigned char*)&plain_text[0]);

    for (size_t i = 0; i < m_blocks; i++){

        nonce_temp[0]=nonce_512; 
        
        AES_encrypt_512(nonce_temp[0], &nonce_temp[0], keys_0, 2);
        
        plain_text[i] =_mm512_xor_si512(plain_text[i],nonce_temp[0]);
        
        AES_encrypt_512(plain_text[i], &plain_text[i], keys_512, 10);

        // imprimiArreglo(16,(unsigned char*)&plain_text[i]);
        // imprimiArreglo(16,(unsigned char*)&plain_text[i]+16);
        // imprimiArreglo(16,(unsigned char*)&plain_text[i]+32);
        // imprimiArreglo(16,(unsigned char*)&plain_text[i]+48);

        S_temp=_mm512_xor_epi64(plain_text[i],S_temp);        
        nonce_512=_mm512_add_epi64(nonce_512, sum_nonce);

    }


    S.bl512=S_temp;
    for (size_t i = 0; i < 4; i++){
        Tag=_mm_xor_si128(Tag,S.bl128[i]);
    }

    AES_encrypt(Tag, &Tag, keys_128, 10);
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