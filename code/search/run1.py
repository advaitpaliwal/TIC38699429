# entry point guarded so TLS spawn workers (macOS) do not re-run the body
import os,sys
if __name__=='__main__':
    sys.argv[0]=os.path.abspath(__file__)
    exec(compile(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'run1_body.py')).read(),'run1_body.py','exec'))
