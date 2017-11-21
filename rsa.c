#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gmp.h>


int a[31];


// Extended Euclidean Algorithm *******************************************
int extended_gcd(mpz_t c, mpz_t a, mpz_t b) {

	mpz_t x, lastx, y, lasty, temp, quotient, temp2;
    
	mpz_init(temp); mpz_init(temp2); mpz_init(quotient); mpz_init_set_ui(x,0); mpz_init_set_ui(lastx,1);mpz_init_set_ui(lasty,0);   
	mpz_init_set_ui(y,1);

    
	while(!(mpz_cmp_ui(b, 0)==0)){
        	mpz_set(temp, b);
	        mpz_tdiv_q(quotient, a, b);
		mpz_mod(b, a, b);
		mpz_set(a, temp);
	       	mpz_set(temp, x);
		mpz_mul(temp2, quotient, x);
		mpz_sub(x, lastx, temp2);
	        mpz_set(lastx, temp);
	        mpz_set(temp, y);
	        mpz_mul(temp2, quotient, y);
		mpz_sub(y, lasty, temp2);
	        mpz_set(lasty, temp);                 
	}
    
	mpz_set(c, lastx); 
	mpz_clear(x); mpz_clear(y); mpz_clear(lastx); mpz_clear(lasty); mpz_clear(temp); mpz_clear(temp2); mpz_clear(quotient);
    
	return 0;		    
}


// Miller - Rabin Primality Test *******************************************
int miller_rabin (mpz_t n){

	int c=0,i=0,flag1=0,flag2=0,flag4=1;

	mpz_t d;
	mpz_t s;
	mpz_t l;
	mpz_t r;
	mpz_t temp1;
	mpz_t temp2;
	mpz_t a_temp;
	mpz_t two;

	mpz_init_set_ui(two,2);
	mpz_init(d);
	mpz_init(s);
	mpz_init(l);
	mpz_init(r);
	mpz_init(temp1);
	mpz_init(a_temp);
	mpz_init(temp2);

	/*a[1]=2; a[2]=3; a[3]=5; a[4]=7; a[5]=11; a[6]=13; a[7]=17; a[8]=19; a[9]=23; a[10]=29;
	a[11]=31;a[12]=37;a[13]=41;a[14]=43;a[15]=47;a[16]=53;a[17]=59;a[18]=61;a[19]=67;a[20]=71;
	a[21]=73;a[22]=79; a[23]=83; a[24]=89; a[25]=97; a[26]=101; a[27]=103; a[28]=107; a[29]=109; a[30]=113;*/

	// Factorization of n to (2^s)*d
	mpz_sub_ui(d, n, 1);

	do {
		mpz_mod_ui(temp1,d,2);
		if(mpz_cmp_ui(temp1,0)==0){
			mpz_tdiv_q_ui(d,d,2);
			mpz_add_ui(s, s, 1);
		} else{
			flag1=1;
		}
	} while(flag1==0);

	//printf("Done first step!\ns = %d and d = %d\n",s,d);

	// End of Factorization

	c=1;

	// loop for different values of a
	do {

		flag2=0;
		mpz_set_ui(a_temp,a[c]);
		mpz_powm(temp1,a_temp,d,n);

		if(mpz_cmp_ui(temp1,1)==0){
			flag2=1;
		}

		mpz_set_ui(r,0);

		if(flag2==0) do {	
			mpz_powm(temp1,two,r,n);
			mpz_mul(l,temp1,d);
			mpz_powm(temp1,a_temp,l,n);
			mpz_sub_ui(temp2,n,1);
	
			if(mpz_cmp(temp1,temp2)==0){
				flag2=1;

			} else{
				mpz_add_ui(r, r, 1);
			}
	
		} while((flag2==0)&&(mpz_cmp(r,s)<0));

		if(flag2==0) flag4=0;
		if(flag4==0) break;
		c++;

	} while(c<31);


	// clearing section
	mpz_clear(d);
	mpz_clear(s);
	mpz_clear(two);
	mpz_clear(l);
	mpz_clear(r);
	mpz_clear(temp1);
	mpz_clear(temp2);
	mpz_clear(a_temp);

	if(flag4==1){
		return 1;
	} else {
		return 0;
	}
}


// key generation function *******************************************
int key_generation(){

	int efound = 0;

	mpz_t p, q, n, e, d, gcd, phi, one, two;
	mp_bitcnt_t bits;
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui (state, time(NULL));

	// initializations
	mpz_init(p);
	mpz_init(q);
	mpz_init(n);
	mpz_init(e);
	mpz_init(d);
	mpz_init(gcd);
	mpz_init(phi);
	mpz_init(two);
	mpz_init(one);

	mpz_set_ui(one,1);
	mpz_set_ui(two,2);
	bits = 512;

	// selecting a random 512 bit prime for p
	do{
		mpz_urandomb(p, state, bits);
	} while (miller_rabin(p) ==0);
	//mpz_out_str(stdout, 10, p);
	//printf("\n");

	// selecting a random 512 bit prime for q
	do{
		mpz_urandomb(q, state, bits);
	} while (miller_rabin(q) ==0);
	//mpz_out_str(stdout, 10, q);
	//printf("\n");

	// computing n
	mpz_mul(n,p,q);
	printf("Public Key = \n");
	mpz_out_str(stdout, 10, n);
	printf("\n\n");

	// computing Φ(n)
	mpz_sub(p, p, one);
	mpz_sub(q, q, one);
	mpz_mul(phi,p,q);

	// choosing a value for e
	mpz_set_ui(e,3);
	mpz_gcd(gcd,e,phi);
	if(mpz_cmp(gcd,one)==0) efound=1;
	if(efound==0){	
		mpz_set_ui(e,17);
	        mpz_gcd(gcd,e,phi);
		if(mpz_cmp(gcd,one)==0) efound=1;
	}
	if(efound ==0){
		mpz_set_ui(e,65537);
        	mpz_gcd(gcd,e,phi);
		if(mpz_cmp(gcd,one)==0) efound=1;
	}
	if (efound==0){
		mpz_set_ui(e,3);
		do { 	mpz_add(e,e,two);
	        	mpz_gcd(gcd,e,phi);
			if(mpz_cmp(gcd,one)==0) efound=1;
		} while (efound ==0);
	}

	printf("e = ");
	mpz_out_str(stdout, 10, e);
	printf("\n\n");

	// computing d
	mpz_t etemp, phitemp;
	mpz_init_set(phitemp, phi);
	mpz_init_set(etemp, e);
	extended_gcd(d,etemp,phitemp);
	mpz_clear(phitemp); mpz_clear(etemp);
	if (mpz_cmp_ui(d,0) < 0) {
		mpz_add(d, d, phi);
	}
	printf("Private Key =\n");
	mpz_out_str(stdout, 10, d);
	printf("\n");

	// deallocations
	mpz_clear(p);
	mpz_clear(q);
	mpz_clear(n);
	mpz_clear(e);
	mpz_clear(d);
	mpz_clear(gcd);
	mpz_clear(phi);
	mpz_clear(one);
	mpz_clear(two);
	//mpz_clear(last);

	//system("PAUSE");
	return 0;
}


// Number Encryption Function *******************************************
int encryption(){

	mpz_t e, n, plain, cipher;
	mpz_init(e); mpz_init(n); mpz_init(plain); mpz_init(cipher);

	printf("Please, enter the public key (n,e).\n");
	mpz_inp_str(n, stdin, 10);
	mpz_inp_str(e, stdin, 10);
	printf("Now enter the number you want to encrypt.\n");
	mpz_inp_str(plain, stdin, 10);
	mpz_powm(cipher, plain, e, n);

	printf("The encrypted message is: ");
	mpz_out_str(stdout, 10, cipher);
	printf("\n");
	
	mpz_clear(e); mpz_clear(n); mpz_clear(plain); mpz_clear(cipher);	
	//system("PAUSE");
	return 1;
}


// Number Decryption Function *******************************************
int decryption(){
	mpz_t n, key, cipher, plain;
	mpz_init(n); mpz_init(key); mpz_init(cipher); mpz_init(plain);
	printf("Please, enter your public and your private key.\n");
	mpz_inp_str(n, stdin, 10);
	mpz_inp_str(key, stdin, 10);
	printf("Now enter the number you want to decrypt.\n");
	mpz_inp_str(cipher, stdin, 10);
	mpz_powm(plain, cipher, key, n);
	printf("The decrypted message is : ");
	mpz_out_str(stdout, 10, plain);
	printf("\n");
	mpz_clear(n); mpz_clear(key); mpz_clear(cipher); mpz_clear(plain);
	//system("PAUSE");
	return 1;
}


// Text Encryption Function *******************************************
int textEncryption(){
	mpz_t n, e, num ,enc;
        mpz_init(e); mpz_init(n); mpz_init(enc);
	mpz_init(num);
		
	char ch;
	int i =0;
	
	srand((unsigned)time(NULL));	
	printf("Please, enter the public key (n,e).\n");
	mpz_inp_str(n, stdin, 10);
	mpz_inp_str(e, stdin, 10);
	printf("Now enter the text message you want to encrypt.\n");

	getchar();

	ch=(char)(rand()%26+100);
	while(ch!='\n'){
	
		mpz_set_ui(num,0);	
		ch=(char)(rand()%26+100);
		while((ch!=' ')&&(ch!='\n')){
		
			mpz_mul_ui(num,num,1000);
			mpz_add_ui(num,num,(int)ch);
			//printf("%c\t%d\n",ch,(int)ch);
			ch=getc(stdin);
	               if (i == 0){
				printf("The encrypted message is : \n");i=1;
			}		
		}

		mpz_powm(enc, num, e, n);
		mpz_out_str(stdout, 10, enc); printf(" ");

	
	}
	printf("\n");
	mpz_clear(num); mpz_clear(n); mpz_clear(e); mpz_clear(enc);
	return 1;
}


// Text Decryption Function *******************************************
int textDecryption() {

	mpz_t n, d, word ,enc,plain,temp,temp2,ten,pow_count,mp;
        mpz_init(d); mpz_init(n); mpz_init(enc);
	mpz_init(word);mpz_init(plain);mpz_init(temp);mpz_init(mp);
	mpz_init(temp2);mpz_init_set_ui(ten,10);mpz_init(pow_count);
	
	char ch,ch2;
	int i =0,flag1,count,num;
	
	srand((unsigned)time(NULL));	
	printf("Please, enter your public and private key.\n");
	mpz_inp_str(n, stdin, 10);
	mpz_inp_str(d, stdin, 10);
	printf("Now enter the text message you want to decrypt.\n");

	ch=' ';
	while(ch==' '){

		mpz_inp_str(word, stdin, 10);
		ch=getchar();	
		mpz_powm(plain, word, d, n);
		count=0;
		flag1=1;
		mpz_set(temp,plain);
		while(flag1==1){
		
			mpz_tdiv_q_ui(temp,temp,10);
			if(mpz_cmp_ui(temp,0)==0){
				flag1=0;
			}
			count++;
		
		}
			
		count=count-3;
		mpz_pow_ui(pow_count,ten,count);
		mpz_tdiv_r(plain,plain,pow_count);

		while(count>0){
			count=count-3;
			mpz_pow_ui(pow_count,ten,count);	
			mpz_tdiv_q(mp,plain,pow_count);
			mpz_tdiv_r(plain,plain,pow_count);		
			ch2=mpz_get_ui(mp);
			printf("%c",ch2);		
		}	
		printf(" ");
	}
	
	printf("\n");
	mpz_clear(word); mpz_clear(n); mpz_clear(d); mpz_clear(enc); mpz_clear(plain);mpz_clear(temp);mpz_clear(temp2);
	mpz_clear(ten);mpz_clear(pow_count);mpz_clear(mp);
	return 1;
}


int main() {

	char input2[80],input;
	int valid = 0, valid2;


	do{
		//menu
		system("clear");
		printf("*********************************************************\n");
		printf("*		  RSA CRYPTOSYSTEM MENU:                *\n");
		printf("*                                                       *\n");
		printf("* -> To generate your public & private keys press g.    *\n");
		printf("*                                                       *\n");
		printf("* -> To encrypt a numeric message press e.              *\n");
		printf("*                                                       *\n");
		printf("* -> To decrypt a numeric message press d.              *\n");
		printf("*                                                       *\n");
		printf("* -> To encrypt a text message press j.                 *\n");
		printf("*                                                       *\n");
		printf("* -> To decrypt a text message press k.                 *\n");
		printf("*                                                       *\n");
		printf("* -> To quit press q.                                   *\n");
		printf("*                                                       *\n");
		printf("*********************************************************\n");

		while (valid == 0){

			scanf("%c",&input);

			if(input == 'g') {
				printf("\n-----Generating Keys------\n\n");
				a[1]=2; a[2]=3; a[3]=5; a[4]=7; a[5]=11; a[6]=13; a[7]=17; a[8]=19; a[9]=23; a[10]=29;
				a[11]=31;a[12]=37;a[13]=41;a[14]=43;a[15]=47;a[16]=53;a[17]=59;a[18]=61;a[19]=67;a[20]=71;
				a[21]=73;a[22]=79; a[23]=83; a[24]=89; a[25]=97; a[26]=101; a[27]=103; a[28]=107; a[29]=109; a[30]=113;

				key_generation(); 
				valid = 1;
	
				printf("\nPress Enter to return to main Menu.");
				getchar();getchar();
	
			} else if (input == 'e') {
				valid = 1; 
				encryption();
				printf("\nPress Enter to return to main Menu.");
				getchar();getchar();
			} else if (input == 'd') {
				valid = 1;
			        decryption();
				printf("\nPress Enter to return to main Menu.");
				getchar();getchar();
			} else if (input == 'j') {
				valid = 1;
				textEncryption();
				printf("\nPress Enter to return to main Menu.");
				getchar();
			} else if (input == 'k') {
				valid = 1;
				textDecryption();
				printf("\nPress Enter to return to main Menu.");
				getchar();
			} else if (input =='q') {
				valid = 1;  
				return 1;
			} else {	
				if(input == '\n') getchar();
				printf("Invalid Input. Please, try again.\n");
				getchar();
			}

		}
	valid=0;

	} while(input !='q');

	return 1;
}


