#!/usr/bin/env python3.6
import sys
from varif import Varif

if __name__ == "__main__":
    argValmap={}
    n=1
    while n <= len(sys.argv[1:]):
        if n == len(sys.argv[1:]):
            if sys.argv[n][0]!="-":
                print("Wrong last argument %s"%sys.argv[n])
                raise(Exception)
            argValmap[sys.argv[n].lstrip("-")]=None
            n+=1
        elif sys.argv[n][0]=="-" and sys.argv[n+1][0]!="-":
            argValmap[sys.argv[n].lstrip("-")]=sys.argv[n+1]
            n+=2
        elif sys.argv[n][0]=="-":
            argValmap[sys.argv[n].lstrip("-")]=None
            n+=1

    Varif(**argValmap)
