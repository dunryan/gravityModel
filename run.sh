if [ -f "satSim.exe" ]; then
    rm "satSim.exe"
    echo "Executable removed."
fi    
g++ *.cpp -o satSim.exe -g
chmod u+x satSim.exe
./satSim.exe

python3 ./plotter_sat.py