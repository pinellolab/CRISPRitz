#include "convert.h"

char convert(char c) {    //fare stessa cosa con la funzione sotto: aggiungere case minuscoli
    switch(c){
        case 'A': c = 0x1; break;
		case 'C': c =  0x2; break;
		case 'G': c =  0x4; break;
		case 'T': c =  0x8; break;
		//case 'N': c =  0x0F; break;
		case '_': c =  0x0F; break;  //inizialmente 0x03, se lo metto al posto di N 0x0F
		case 'R': c =  0x05; break;
		case 'Y': c =  0x0A; break;
		case 'S': c =  0x06; break;
		case 'W': c =  0x09; break;
		case 'K': c =  0x0C; break;
		case 'M': c =  0x03; break; //attenzione, è uguale a '_',  cambio '_' con N
		case 'B': c =  0x0E; break;
		case 'D': c =  0x0D; break;
		case 'H': c =  0x0B; break;
		case 'V': c =  0x07; break;
		default : std::cout << "Golden wind";c =  0x0; break;
    }
    return c;
}


std::string convertForBitset(char c){
	std::string result;
	switch(c){
		case 'A': case 'a': result = "0001"; break;
		case 'C': case 'c': result = "0010"; break;
		case 'G': case 'g': result = "0100"; break;
		case 'T': case 't': result = "1000"; break;
		//case 'N': c =  0x0F; break;
		case '_': result = "1111"; break;  //inizialmente 0x03, se lo metto al posto di N 0x0F
		case 'R': case 'r': result = "0101"; break;
		case 'Y': case 'y': result = "1010"; break;
		case 'S': case 's': result = "0110"; break;
		case 'W': case 'w': result = "1001"; break;
		case 'K': case 'k': result = "1100"; break;
		case 'M': case 'm': result = "0011"; break; //attenzione, è uguale a '_',  cambio '_' con N
		case 'B': case 'b': result = "1110"; break;
		case 'D': case 'd': result = "1101"; break;
		case 'H': case 'h': result = "1011"; break;
		case 'V': case 'v': result = "0111"; break;
		default : result = "1111"; break; //indico la N, attenzione forse mi serve un altro per '-'
	}
	return result;
}

char convert(std::bitset<4> nuc_bit){
	auto nuc_int = nuc_bit.to_ulong();
	switch(nuc_int) {
		case 1: return 'A';
		case 2: return 'C';
		case 4: return 'G';
		case 8: return 'T';
		case 5: return 'R';
		case 10: return 'Y';
		case 6: return 'S';
		case 9: return 'W';
		case 12: return 'K';
		case 3: return 'M';
		case 14: return 'B';
		case 13: return 'D';
		case 11: return 'H';
		case 7: return 'V';
		case 0: return '-';
		default: //std::cout << "Il numero: " << nuc_int << " "; break;
		break;
	}
	return 'N';
}