/*
 Functions used by wsprsim
 */

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>

#include "wsprsim_utils.h"
#include "wsprd_utils.h"
#include "nhash.h"
#include "fano.h"


char get_locator_character_code(char ch) {
    if( ch >=48 && ch <=57 ) { //0-9
        return ch-48;
    }
    if( ch == 32 ) {  //space
        return 36;
    }
    if( ch >= 65 && ch <= 82 ) { //A-Z
        return ch-65;
    }
    return -1;
}

char get_callsign_character_code(char ch) {
    if( ch >=48 && ch <=57 ) { //0-9
        return ch-48;
    }
    if( ch == 32 ) {  //space
        return 36;
    }
    if( ch >= 65 && ch <= 90 ) { //A-Z
        return ch-55;
    }
    return -1;
}

long unsigned int pack_grid4_power(char *grid4, int power) {
    long unsigned int m;

    m=(179-10*grid4[0]-grid4[2])*180+10*grid4[1]+grid4[3];
    m=m*128+power+64;
    return m;
}

long unsigned int pack_call(char *callsign) {
    int i;
    long unsigned int n;
    char call6[6];
    memset(call6,32,sizeof(char)*6);
    // callsign is 6 characters in length. Exactly.
    int call_len = strlen(callsign);
    if( call_len > 6 ) {
        return 0;
    }
    if( isdigit(*(callsign+2)) ) {
        for (i=0; i<6; i++) {
            if( callsign[i] == 0 ) {
                call6[i]=32;
            } else {
                call6[i]=*(callsign+i);
            }
        }
    } else if( isdigit(*(callsign+1)) ) {
        for (i=0; i<6; i++) {
            if( i==0 || callsign[i-1]==0 ) {
                call6[i]=32;
            } else {
                call6[i]=*(callsign+i-1);
            }
        }
    }
    for (i=0; i<6; i++) {
        call6[i]=get_callsign_character_code(call6[i]);
    }
    n = call6[0];
    n = n*36+call6[1];
    n = n*10+call6[2];
    n = n*27+call6[3]-10;
    n = n*27+call6[4]-10;
    n = n*27+call6[5]-10;
    return n;
}

void pack_prefix(char *callsign, int32_t *n, int32_t *m, int32_t *nadd ) {
    int i;
    char *call6;
    call6=malloc(sizeof(char)*6);
    memset(call6,32,sizeof(char)*6);
    int i1=strcspn(callsign,"/");

    if( callsign[i1+2] == 0 ) {
        //single char suffix
        for (i=0; i<i1; i++) {
            call6[i]=callsign[i];
        }
        *n=pack_call(call6);
        *nadd=1;
        int nc = callsign[i1+1];
        if( nc >= 48 && nc <= 57 ) {
            *m=nc-48;
        } else if ( nc >= 65 && nc <= 90 ) {
            *m=nc-65+10;
        } else {
            *m=38;
        }
        *m=60000-32768+*m;
    } else if( callsign[i1+3]==0 ) {
        //two char suffix
        for (i=0; i<i1; i++) {
            call6[i]=callsign[i];
        }
        *n=pack_call(call6);
        *nadd=1;
        *m=10*(callsign[i1+1]-48)+(callsign[i1+2]-48);
        *m=60000 + 26 + *m;
    } else {
        char* pfx=strtok(callsign,"/");
        call6=strtok(NULL," ");
        *n=pack_call(call6);
        int plen=strlen(pfx);
        if( plen ==1 ) {
            *m=36;
            *m=37*(*m)+36;
        } else if( plen == 2 ) {
            *m=36;
        } else {
            *m=0;
        }
        for (i=0; i<plen; i++) {
            int nc = callsign[i];
            if( nc >= 48 && nc <= 57 ) {
                nc=nc-48;
            } else if ( nc >= 65 && nc <= 90 ) {
                nc=nc-65+10;
            } else {
                nc=36;
            }
            *m=37*(*m)+nc;
        }
        *nadd=0;
        if( *m > 32768 ) {
            *m=*m-32768;
            *nadd=1;
        }
    }
}

void interleave(unsigned char *sym) {
    unsigned char tmp[162];
    unsigned char p, i, j;

    p=0;
    i=0;
    while (p<162) {
        j=((i * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
        if (j < 162 ) {
            tmp[j]=sym[p];
            p=p+1;
        }
        i=i+1;
    }
    for (i=0; i<162; i++) {
        sym[i]=tmp[i];
    }
}

int get_wspr_channel_symbols(char* rawmessage, char* hashtab, unsigned char* symbols) {
    int m=0, n=0, ntype=0;
    int i, j, ihash;
    unsigned char pr3[162]= {
        1,1,0,0,0,0,0,0,1,0,0,0,1,1,1,0,0,0,1,0,
        0,1,0,1,1,1,1,0,0,0,0,0,0,0,1,0,0,1,0,1,
        0,0,0,0,0,0,1,0,1,1,0,0,1,1,0,1,0,0,0,1,
        1,0,1,0,0,0,0,1,1,0,1,0,1,0,1,0,1,0,0,1,
        0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,0,0,0,1,0,
        0,0,0,0,1,0,0,1,0,0,1,1,1,0,1,1,0,0,1,1,
        0,1,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,1,1,
        0,0,0,0,0,0,0,1,1,0,1,0,1,1,0,0,0,1,1,0,
        0,0
    };
    int nu[10]= {0,-1,1,0,-1,2,1,0,-1,1};
    char *callsign, *grid, *powstr;
    char grid4[5], message[23];

    memset(message,0,sizeof(char)*23);
    i=0;
    while ( rawmessage[i] != 0 && i<23 ) {
        message[i]=rawmessage[i];
        i++;
    }

    int i1=strcspn(message," ");
    int i2=strcspn(message,"/");
    int i3=strcspn(message,"<");
    int i4=strcspn(message,">");
    int mlen=strlen(message);

    // Use the presence and/or absence of "<" and "/" to decide what
    // type of message. No sanity checks! Beware!

    if( (i1>3) & (i1<7) & (i2==mlen) & (i3==mlen) ) {
        // Type 1 message: K9AN EN50 33
        //                 xxnxxxx xxnn nn
        callsign = strtok(message," ");
        grid = strtok(NULL," ");
        powstr = strtok(NULL," ");
        int power = atoi(powstr);
        n = pack_call(callsign);

        for (i=0; i<4; i++) {
            grid4[i]=get_locator_character_code(*(grid+i));
        }
        m = pack_grid4_power(grid4,power);

    } else if ( i3 == 0 && i4 < mlen ) {
        // Type 3:      <K1ABC> EN50WC 33
        //          <PJ4/K1ABC> FK52UD 37
        // send hash instead of callsign to make room for 6 char grid.
        // if 4-digit locator is specified, 2 spaces are added to the end.
        callsign=strtok(message,"<> ");
        grid=strtok(NULL," ");
        powstr=strtok(NULL," ");
        int power = atoi(powstr);
        if( power < 0 ) power=0;
        if( power > 60 ) power=60;
        power=power+nu[power%10];
        ntype=-(power+1);
        ihash=nhash(callsign,strlen(callsign),(uint32_t)146);
        m=128*ihash + ntype + 64;

        char grid6[6];
        memset(grid6,32,sizeof(char)*6);
        j=strlen(grid);
        for(i=0; i<j-1; i++) {
            grid6[i]=grid[i+1];
        }
        grid6[5]=grid[0];
        n=pack_call(grid6);
    } else if ( i2 < mlen ) {  // just looks for a right slash
        // Type 2: PJ4/K1ABC 37
        callsign=strtok(message," ");
        if( strlen(callsign) < i2 ) return 0; //guards against pathological case
        powstr=strtok(NULL," ");
        int power = atoi(powstr);
        if( power < 0 ) power=0;
        if( power > 60 ) power=60;
        power=power+nu[power%10];
        int n1, ng, nadd;
        pack_prefix(callsign, &n1, &ng, &nadd);
        ntype=power + 1 + nadd;
        m=128*ng+ntype+64;
        n=n1;
    } else {
        return 0;
    }

    // pack 50 bits + 31 (0) tail bits into 11 bytes
    unsigned char it, data[11];
    memset(data,0,sizeof(char)*11);
    it=0xFF & (n>>20);
    data[0]=it;
    it=0xFF & (n>>12);
    data[1]=it;
    it=0xFF & (n>>4);
    data[2]=it;
    it= ((n&(0x0F))<<4) + ((m>>18)&(0x0F));
    data[3]=it;
    it=0xFF & (m>>10);
    data[4]=it;
    it=0xFF & (m>>2);
    data[5]=it;
    it=(m & 0x03)<<6 ;
    data[6]=it;
    data[7]=0;
    data[8]=0;
    data[9]=0;
    data[10]=0;


    // make sure that the 11-byte data vector is unpackable
    // unpack it with the routine that the decoder will use and display
    // the result. let the operator decide whether it worked.
//    char hashtab[32768][13];
//    memset(hashtab,0,sizeof(char)*32768*13);

    char *check_call_loc_pow, *check_callsign, *call, *loc, *pwr;
    check_call_loc_pow=malloc(sizeof(char)*23);
    check_callsign=malloc(sizeof(char)*13);
    call=malloc(sizeof(char)*13);
    loc=malloc(sizeof(char)*7);
    pwr=malloc(sizeof(char)*3);
    signed char check_data[11];
    memcpy(check_data,data,sizeof(char)*11);
    unpk_(check_data,hashtab,check_call_loc_pow,call,loc,pwr,check_callsign);
//    printf("Will decode as: %s\n",check_call_loc_pow);

    unsigned int nbytes=11; // The message with tail is packed into 11 bytes.
    unsigned int nencoded=162;
    unsigned char channelbits[nencoded];
    memset(channelbits,0,sizeof(char)*nencoded);

    encode(channelbits,data,nbytes);

    interleave(channelbits);

    for (i=0; i<162; i++) {
        symbols[i]=2*channelbits[i]+pr3[i];
    }

    return 1;
}
