from subprocess import Popen, PIPE, STDOUT
import sys
import tarfile

def run_shell_command(command):
    with Popen(
        command, shell=True, universal_newlines=True, bufsize=1, stdout=PIPE, stderr=STDOUT
    ) as p:
        for out in p.stdout:
            print(out, end='')
    if p.returncode != 0:
        print('Exit code: %d' % p.returncode)
        sys.exit(1)


def untar(tar_file, untar_to):
    print('untar files to %s ...' % untar_to, flush=True, end='')
    # extract files from tar ball
    try:
        with tarfile.open(tar_file, 'r:gz') as tar:
            tar.extractall(path=untar_to)
        print('Done.')
    except Exception as e:
        print('/nError when extracting files: %s' % str(e))
        sys.exit(1)

def run_templates_in_shell(list_of_templates, mapping):
    for template in list_of_templates:
        run_shell_command(template.substitute(mapping))
