from subprocess import Popen, PIPE, STDOUT
import sys
import os
import tarfile


def run_shell_command(command):
    print('+' + command, flush=True)
    with Popen(
        command, shell=True, universal_newlines=True, bufsize=1, stdout=PIPE, stderr=STDOUT
    ) as p:
        for out in p.stdout:
            print(out, end='')
    if p.returncode != 0:
        sys.exit('Exit code: %d' % p.returncode)


def untar(tar_file, untar_to):
    print('untar files to %s ...' % untar_to, flush=True, end='')
    # extract files from tar ball
    try:
        with tarfile.open(tar_file, 'r:gz') as tar:
            tar.extractall(path=untar_to)
        print('Done.')
    except Exception as e:
        sys.exit('Error: unexpected exception. When extracting files: %s' % str(e))


def run_templates_in_shell(list_of_templates, mapping):
    for template in list_of_templates:
        run_shell_command(template.substitute(mapping))


def mkdir(dir_path):
    if not os.path.isdir(dir_path):
        try:
            os.makedirs(dir_path)
        except Exception as e:
            sys.exit('Error: unexpected exception. failed to create the directory %s: %s' % (dir_path, str(e)))
