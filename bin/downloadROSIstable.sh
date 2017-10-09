# One click install script dumux-rosi

# dumux-multidimension
echo "$(tput setaf 2)Downloading dumux-multidimension...$(tput setaf 7)"
git clone -b feature/embeddedcoupling-Reuse-source-sink-value https://git.iws.uni-stuttgart.de/dumux-appl/dumux-multidimension.git

# the core modules
echo "$(tput setaf 2)Downloading DUNE core modules...$(tput setaf 7)"
for MOD in common geometry grid localfunctions istl; do
    if [ ! -d "dune-$MOD" ]; then
        git clone -b releases/2.5 https://gitlab.dune-project.org/core/dune-$MOD.git
    else
        echo "Skip cloning dune-$MOD because the folder already exists."
        cd dune-$MOD
        git checkout releases/2.5
        cd ..
    fi
done

git clone -b releases/2.5 https://gitlab.dune-project.org/staging/dune-uggrid.git

git clone -b releases/2.5 https://gitlab.dune-project.org/extensions/dune-foamgrid.git


echo "$(tput setaf 2)Downloading DUMUX ...$(tput setaf 7)"
git clone -b releases/2.11 https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git


# dumux-rosi
echo "$(tput setaf 2)Downloading dumux-rosi ...$(tput setaf 7)"
git clone https://github.com/Plant-Root-Soil-Interactions-Modelling/dumux-rosi.git




