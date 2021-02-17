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


std::string convertCharToBitset(char c){
	std::string result;
	switch(c){
		case 'A': case 'a': result = "0001"; break;
		case 'C': case 'c': result = "0010"; break;
		case 'G': case 'g': result = "0100"; break;
		case 'T': case 't': result = "1000"; break;
		case 'R': case 'r': result = "0101"; break;
		case 'Y': case 'y': result = "1010"; break;
		case 'S': case 's': result = "0110"; break;
		case 'W': case 'w': result = "1001"; break;
		case 'K': case 'k': result = "1100"; break;
		case 'M': case 'm': result = "0011"; break; 
		case 'B': case 'b': result = "1110"; break;
		case 'D': case 'd': result = "1101"; break;
		case 'H': case 'h': result = "1011"; break;
		case 'V': case 'v': result = "0111"; break;
		case 'N': case 'n': result = "1111"; break;
		case '-': 			result = "0000"; break;
		default : std::cerr << "The character (" << c << ") is not part of the IUPAC nucleotide nomenclature" << std::endl; 
				  result = "1111"; break;
	}
	return result;
}

char convertBitsetToChar(std::bitset<4> nuc_bit){

	if(nuc_bit == std::bitset<4>("0001"))
		return 'A';
	else if(nuc_bit == std::bitset<4>("0010"))
		return 'C';
	else if(nuc_bit == std::bitset<4>("0100"))
		return 'G';
	else if(nuc_bit == std::bitset<4>("1000"))
		return 'T';
	else if(nuc_bit == std::bitset<4>("0101"))
		return 'R';
	else if(nuc_bit == std::bitset<4>("1010"))
		return 'Y';
	else if(nuc_bit == std::bitset<4>("0110"))
		return 'S';
	else if(nuc_bit == std::bitset<4>("1001"))
		return 'W';
	else if(nuc_bit == std::bitset<4>("1100"))
		return 'K';
	else if(nuc_bit == std::bitset<4>("0011"))
		return 'M';
	else if(nuc_bit == std::bitset<4>("1110"))
		return 'B';
	else if(nuc_bit == std::bitset<4>("1101"))
		return 'D';
	else if(nuc_bit == std::bitset<4>("1011"))
		return 'H';
	else if(nuc_bit == std::bitset<4>("0111"))
		return 'V';
	else if(nuc_bit == std::bitset<4>("0000"))
		return '-';
	else if(nuc_bit == std::bitset<4>("1111"))
		return 'N';
	else
		std::cout << "il nucbit inserito è: " << nuc_bit << std::endl;
	return '?';
}