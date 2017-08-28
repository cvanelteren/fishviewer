from pip import main


def packageInstaller(packages):
    from pip import main

    for p in packages:
            main(['install', p])

if __name__ == '__main__':
    packages = 'PyQt5 tqdm scipy sklearn h5py numpy matplotlib spyder seaborn'

    packageInstaller(packages.split(' '))
