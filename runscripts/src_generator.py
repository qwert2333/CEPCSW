#!/usr/bin/env python

import os
from time import time
from string import Template


def generate():

    template_file = open(
        r'/afs/ihep.ac.cn/users/h/huangliheng/work/CEPCSW/runscripts/read.py', 'r')
    tmpl = Template(template_file.read())

    if not os.path.exists('/afs/ihep.ac.cn/users/h/huangliheng/work/CEPCSW/runscripts/bashes/'):
        os.makedirs(
            '/afs/ihep.ac.cn/users/h/huangliheng/work/CEPCSW/runscripts/bashes/')

    filenames = []
    for index in range(10):

        num = index
        rndm = int(time()*10000)

        lines = []
        lines.append(tmpl.substitute(NUM=num))

        output_file_name = 'read_%d.py' % (rndm)
        filenames.append(output_file_name)
        output_file = open(
            '/afs/ihep.ac.cn/users/h/huangliheng/work/CEPCSW/runscripts/bashes/' + output_file_name, 'w')

        output_file.writelines(lines)
        output_file.close()

        print('generate job file : ', output_file_name)

        with open('/afs/ihep.ac.cn/users/h/huangliheng/work/CEPCSW/runscripts/bashes/filenames.txt', 'w') as f:
            for filename in filenames:
                f.write("%s\n" % filename)


if __name__ == '__main__':

    generate()
