

// icc -o extractnonfunc extractnonfunc.c

#include <stdio.h>
#include <stdlib.h>



// just use stdin stdout
int main(void)
{
  int countbraceopen,countparopen;
  int commentopen; // /* ... */
  char ch,ch2;
  char lastch='A';


  countbraceopen=countparopen=commentopen=0;

  while(!feof(stdin)){

    lastch=ch;
    ch=fgetc(stdin);

    if((ch=='\n' || ch=='\r')){
      // see if next is open brace corresponding to new function
      ch2=fgetc(stdin);
      if(ch2=='{'){
	if(countbraceopen==0 && countparopen==0) fprintf(stdout," FUNCBEGIN ");	
	countbraceopen++;
      }
      else{
	// then pop back one car and copy output
	//	fputc(ch,stdout);
	//	fseek(stdin,-1,SEEK_CUR);

	if(ch2=='{'){
	  if(countbraceopen==0) fprintf(stdout," BRACEBEGIN ");
	  countbraceopen++;
	}
	else if(ch2=='}') countbraceopen--;
	else if(ch2=='('){
	  //      fputc(ch2,stdout);
	  if(countbraceopen==0 && countparopen==0) fprintf(stdout," PARBEGIN ");
	  countparopen++;
	}
	else if(ch2==')'){
	  //      fputc(ch2,stdout);
	  countparopen--;
	}
	else if(ch=='/' && ch2=='*'){
	  commentopen++;
	}
	else if(ch=='*' && ch2=='/'){
	  commentopen--;
	}
	else{
	  if(countbraceopen==0 && countparopen==0 && commentopen==0){
	    fputc(ch2,stdout);
	    //	if(countparopen==0 || (countparopen==1 && (ch2!='\n' || ch2!='\r'))) fputc(ch2,stdout);
	  }
	}

      }
    }
    else if(ch=='{'){
      if(countbraceopen==0) fprintf(stdout," BRACEBEGIN ");
      countbraceopen++;
    }
    else if(ch=='}') countbraceopen--;
    else if(ch=='('){
      //      fputc(ch,stdout);
      if(countbraceopen==0 && countparopen==0) fprintf(stdout," PARBEGIN ");
      countparopen++;
    }
    else if(ch==')'){
      //      fputc(ch,stdout);
      countparopen--;
    }
    else if(lastch=='/' && ch=='*'){
      commentopen++;
    }
    else if(lastch=='*' && ch=='/'){
      commentopen--;
    }
    else{
      if(countbraceopen==0 && countparopen==0 && commentopen==0){
	fputc(ch,stdout);
	//	if(countparopen==0 || (countparopen==1 && (ch!='\n' || ch!='\r'))) fputc(ch,stdout);
      }
    }
  }// end over file


  return(0);
}


// then use: ./extractnonfunc < coord.c | grep -v "^void " | grep -v "^//"| grep -v "^/\*" | grep -v "^#" | grep -v "FUNCBEGIN" | less
// See codenonfunc.sh

