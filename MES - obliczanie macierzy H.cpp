// MES - obliczanie macierzy H.cpp : Ten plik zawiera funkcję „main”. W nim rozpoczyna się i kończy wykonywanie programu.
//

#include <iostream>
#include <cmath>

struct Node {
    int id;
    double x, y;
    int warunekBrzegowy = 0;

    /*
    void uzupelnijNode(int identyfikator, double var1, double var2, double h, double b) {
        id = identyfikator;
        x = var1;
        y = var2;
        
        warunekBrzegowy = ustalenieWarunkuBrzegowego(x,y,h,b);        
        //this->warunekBrzegowy = ustalenieWarunkuBrzegowego(x, y, h, b); 
    }
    */

    Node(){
        this->id = 0;
        this->x = 0;
        this->y = 0;
        this->warunekBrzegowy = 0;
    }

    Node(int identyfikator, double var1, double var2, double h, double b) {
        id = identyfikator;
        x = var1;
        y = var2;

        warunekBrzegowy = ustalenieWarunkuBrzegowego(x, y, h, b);
    }


    void wyswietlNode() {
        std::cout << "[" << x << "," << y << "]\t"<<"BC: "<<warunekBrzegowy<<"\t";
    }

    void wyswietlNodeNumer() {
        std::cout << "Node nr: " << id << "[" << x << "," << y << "]\n";
    }

    int ustalenieWarunkuBrzegowego(double var1, double var2, double h, double b) {
        if ((var1 == 0) && (var2 != h)) {
            //std::cout << "x = " << var1 << "\ty = " << var2 << "\n";
            //std::cout << "h = " << h << "\tb = " << b << "\n";
            return 1;
        }
        else if ((var1 != b) && (var2 == 0)) {
            //std::cout << "x = " << var1 << "\ty = " << var2 << "\n";
            //std::cout << "h = " << h << "\tb = " << b << "\n";
            return 1;
        }

        else if ((var1 != b) && (var2 == h)) {
            //std::cout << "x = " << var1 << "\ty = " << var2 << "\n";
            //std::cout << "h = " << h << "\tb = " << b << "\n";
            return 1;
        }
        else if ((var1 == b) && (var2 != h)) {
            //std::cout << "x = " << var1 << "\ty = " << var2 << "\n";
            //std::cout << "h = " << h << "\tb = " << b << "\n";
            return 1;
        }
        else if ((var1 == b) && (var2 == h)) {
            return 1;
        }
        else {
            return 0;
        }
    }

    /*
    void show()
{
    cout << "id: " << id << endl;
    cout << "bc: " << bc << endl;
    cout << "( " << x << ", " << y << " )\n";
}
    */
};

struct Element {
    int value;
    int nN;
    //int NODE[4];

    double** MacierzH;
    double** MacierzHbc;
    double** MacierzC;
    double* wektorP;
    Node* nodes;

    /*
    void uzupelnijElement(int wartosc, int b, int c) {
        value = wartosc;
        NODE[0] = b;
        NODE[1] = b + c;
        NODE[2] = NODE[1] + 1;
        NODE[3] = b + 1;

        MacierzH = new double* [4];
        MacierzHbc = new double* [4];
        MacierzC = new double* [4];
        wektorP = new double[4];

        //Ustawienie wartości początkowych dla macierzy
        for (int i = 0; i < 4; i++) {
            MacierzH[i] = new double[4];
            MacierzHbc[i] = new double[4];
            MacierzC[i] = new double[4];
            wektorP[i] = 0;
            for (int j = 0; j < 4; j++) {
                MacierzH[i][j] = 0;
                MacierzHbc[i][j] = 0;
                MacierzC[i][j] = 0;
            }
        }
    }*/
    
    Element() {}
    
    Element(int wartosc, int b, int c, Node* arr, int nn) {
        value = wartosc;
        nN = nn;
        //NODE[0] = arr[b];
        //NODE[1] = arr[b + c];
        //NODE[2] = arr[NODE[1] + 1];
        //NODE[3] = arr[b + 1];
        nodes = new Node[4];
        nodes[0] = arr[b - 1];
        nodes[1] = arr[b + c - 1];
        nodes[2] = arr[nodes[1].id];
        nodes[3] = arr[b];

        MacierzH = new double* [4];
        MacierzHbc = new double* [4];
        MacierzC = new double* [4];
        wektorP = new double[4];
        
        //Ustawienie wartości początkowych dla macierzy
        for (int i = 0; i < 4; i++) {
            MacierzH[i] = new double[4];
            MacierzHbc[i] = new double[4];
            MacierzC[i] = new double[4];
            wektorP[i] = 0;
            for (int j = 0; j < 4; j++) {
                MacierzH[i][j] = 0;
                MacierzHbc[i][j] = 0;
                MacierzC[i][j] = 0;
            }
        }

        /*
        std::cout << "MacierzHbc: \n";
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                std::cout << MacierzHbc[i][j] << "\t";
            }
            std::cout << "\n";
        }
        */
    }
    
    /*
    void wyswietlElement() {
        std::cout << "Element nr. " << value << " = ";
        for (int i = 0; i < 4; i++) {
            std::cout << NODE[i] << "\t";
        }
        std::cout << "\n";
    }
    */

    
    void wyswietlElement() {
        std::cout << "Element: " << value << "\n";
        std::cout << "Wezly: " << nodes[0].id << ", " << nodes[1].id << ", " << nodes[2].id << ", " << nodes[3].id << "\n";
    }

    /*void show()
    {
        std::cout << "value: " << value << "\n";
        std::cout << "nodes: " << nodes[0].id << ", " << nodes[1].id << ", " << nodes[2].id << ", " << nodes[3].id << "\n";
    }*/
    /*
    
    void showMore()
    {
        show();
        nodes[0].show();
        cout << "\n";
        nodes[1].show();
        cout << "\n";
        nodes[2].show();
        cout << "\n";
        nodes[3].show();
        cout << "\n";
    }


   */

    void wyświetlH() {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                std::cout << MacierzH[i][j] << "\t";
            }                
            std::cout << "\n";
        }
        std::cout << "\n";
    }
    
    void wyświetlHBC() {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                std::cout << MacierzHbc[i][j] << "\t";
            }                
            std::cout << "\n";
        }
        std::cout << "\n";
    }
     

    void obliczMacierzC(int i, double** n, double c, double ro, double** J, int punktyCalkowania, double waga) {
        double wyznacznik = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        std::cout << "Wyznacznik dla obliczania macierzy C: " << wyznacznik << "\n";
        if (punktyCalkowania == 2) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    MacierzC[j][k] += c * ro * wyznacznik * (n[i][j] * n[i][k]);
                    //std::cout << MacierzC[j][k] << "\n";
                }
            }
        }
        else if (punktyCalkowania == 3) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    MacierzC[j][k] += c * ro * wyznacznik * (n[i][j] * n[i][k]) * waga;
                    //std::cout << MacierzC[j][k] << "\n";
                }
            }
        }
    }

    double długość(int surfaceId) {
        int x1 = 0, x2 = 0;
        switch (surfaceId) {
            
            case 0: 
            {//dół
                x1 = 0;
                x2 = 1;
                break;
            }
        
            case 1:
            {//prawa
                x1 = 1;
                x2 = 2;
                break;
            }
        
            case 2:
            {//góra
                x1 = 2;
                x2 = 3;
                break;
            }
        
            case 3:
            {//lewa
                x1 = 3;
                x2 = 0;
                break;
            }
        }

        double xDiff = nodes[x1].x - nodes[x2].x;
        double yDiff = nodes[x1].y - nodes[x2].y;
        //return 0.00015625;
        return sqrt(xDiff * xDiff + yDiff * yDiff);
    }

    void agregacja(double** tab)
    {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                MacierzHbc[i][j] += MacierzH[i][j];
            }                
        }
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                int x = nodes[i].id - 1;
                int y = nodes[j].id - 1;
                tab[y][x] = MacierzHbc[i][j];
            }
        }
    }

    void agregacjaP(double* vec) {
        for (int i = 0; i < 4; i++) {
            vec[nodes[i].id - 1] = wektorP[i];
        }            
    }

    void obliczanieP(double* n1, double* n2, double* n3, double t, int surfaceId, int pc) {
        double det = długość(surfaceId) / 2;
        if (pc == 2) {
            double w1 = 1, w2 = 1;
            for (int i = 0; i < 4; i++) {
                wektorP[i] += 300.0 * ((w1 * n1[i]) * t + (w2 * n2[i]) * t) * det;
                //std::cout << "Test\t wektorP[" << i << "] = " << wektorP[i] << "\n";
            }
        }
        else if (pc == 3) {
            double w1 = 5.0 / 9.0, w2 = 8.0 / 9.0, w3 = w1;
            for (int i = 0; i < 4; i++) {
                wektorP[i] += 300.0 * ((w1 * n1[i]) * t + (w2 * n2[i]) * t + (w1 * n3[i]) * t) * det;
                //std::cout << "Test\t wektorP[" << i << "] = " << wektorP[i] << "\n";
            }
        }
    }

    void agregacjaC(double** tab) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                int x = nodes[i].id - 1;
                int y = nodes[j].id - 1;
                tab[y][x] = MacierzC[i][j];
            }
        }
    }
};

struct Grid {
    double h, b;
    int nH, nB, nN, nE;
    int punktyCalkowania, initialTemperature, simulationStepTime, ambientTemperature, alfa, specificHeat, conductivity, density;
    Node* nodes;
    Element* elements;
    double** nc;
    double** globalneH;
    double** globalneC;
    double* globalnyP;
    double* T0;
    double* koncowyP;

    /*
    Grid(double H, double B, int NH, int NB, int pc) {
        punktyCalkowania = pc;
        h = H;
        b = B;
        nH = NH;
        nB = NB;
        nN = nH * nB;
        nE = (nH - 1) * (nB - 1);
        nodes = new Node[nN];
        elements = new Element[nE];

       
        double dx = b / (nB - 1), dy = h / (nH - 1);
        int k = 0;
        for (int i = 0; i < nB; ++i)
            for (int j = 0; j < nH; ++j)
            {
                nodes[k] = Node(k + 1, dx * i, dy * j, h, b);
                k++;
            }

        int ID = 1;
        for (int i = 0; i < nE; ++i)
        {
            if (i != 0)
                if (i % (nB - 1) == 0)
                    ID++;
            elements[i] = Element(i + 1, ID, nH, nodes, nN);
            ID = elements[i].nodes[3].id;
        }

        globalnyP = new double[nN];
        koncowyP = new double[nN];
        globalneH = new double* [nN];
        globalneC = new double* [nN];
        T0 = new double[nN];
        for (int i = 0; i < nN; ++i)
        {
            globalnyP[i] = 0.0;
            koncowyP[i] = 0.0;
            T0[i] = 100.0;
            globalneH[i] = new double[nN];
            globalneC[i] = new double[nN];
            for (int j = 0; j < nN; ++j)
            {
                globalneH[i][j] = 0.0;
                globalneC[i][j] = 0.0;
            }
        }
    }
    */

    Grid(double wys, double szer, int nodePion, int nodePoziom, int integrationPoints, int tPocz, int krokCzasowy, int tOtocz, int wspAlfa, int cieploWlasciwe, int cond, int gestosc) {
        //przypisanie wartosci
        h = wys;
        b = szer;
        nH = nodePion;
        nB = nodePoziom;
        nN = nH * nB;
        nE = (nH - 1) * (nB - 1);

        punktyCalkowania = integrationPoints;
        initialTemperature = tPocz;
        simulationStepTime = krokCzasowy;
        ambientTemperature = tOtocz;
        alfa = wspAlfa;
        specificHeat = cieploWlasciwe;
        conductivity = cond;
        density = gestosc;

        nodes = new Node[nN];
        elements = new Element[nE];
        //stworzenie node'ow        
        double dx = b / (nodePoziom - 1);
        double dy = h / (nodePion - 1);

        int g = 0;
        for (int i = 0; i < nB; i++) {
            for (int j = 0; j < nH; j++) {
                //nodes[g].uzupelnijNode(g + 1, dx * i, dy * j, h, b);      
                nodes[g] = Node(g + 1, dx * i, dy * j, h, b);
                g++;                
            }
        }


        //stworzenie elementow
        
        int ID = 1;
        for (int i = 0; i < nE; i++) {
            if (i != 0) {
                if (i % (nH - 1) == 0) {
                    ID++;
                }
            }            
            //elements[i].uzupelnijElement(i + 1, ID, nH);
            elements[i] = Element(i + 1, ID, nH, nodes, nN);
            ID = elements[i].nodes[3].id;
        }

        //Tworzenie globalnych wektorów i macierzy
        globalneH = new double*[nN];
        globalneC = new double*[nN];

        for (int i = 0; i < nN; i++) {
            globalneC[i] = new double[nN];
            globalneH[i] = new double[nN];
        }

        globalnyP = new double[nN];
        T0 = new double[nN];
        koncowyP = new double[nN];

        //Ustawienie wartości początkowych dla stworzonych wektorów i macierzy
        for (int i = 0; i < nN; i++) {
            globalnyP[i] = 0.0;
            T0[i] = tPocz;
            koncowyP[i] = 0.0;

            for (int j = 0; j < nN; j++) {
                globalneC[i][j] = 0.0;
                globalneH[i][j] = 0.0;
            }
        }

    }

    void wyswietlElementy() {
        std::cout << "\nElementy:\n";
        for (int i = 0; i < nE; i++) {
            elements[i].wyswietlElement();
        }
        std::cout << "\n";
    }

    void wyswietlNody() {
        std::cout << "\nNode'y\n";
        for (int i = 0; i < nN; i++) {
            nodes[i].wyswietlNodeNumer();
        }
        std::cout << "\n";
    }

    void wyswietlSiatke() {
        std::cout << "\nSiatka:\n";
        for (int i = nH; i > 0; i--) {
            int k = i - 1;
            for (int j = 0; j < nB; j++) {
                nodes[k].wyswietlNode();
                k += nH;
            }
            std::cout << "\n";
        }
        std::cout << "\n";

        
        std::cout<< "Elementy: \n";
        for(int i = 0; i < (nB - 1)*(nH-1); i++){
            elements[i].wyswietlElement();
        }
        std::cout << "\n";

        std::cout << "Node'y:\n";
        for (int i = 0; i < nN; i++) {
            nodes[i].wyswietlNodeNumer();
        }
        std::cout << "\n";
        
    }

    void funkcjaJakobian(double** pochodnaKsi, double** pochodnaEta) {

        double** Jakobian = new double* [2];
        double** JakobianOdwrotny = new double* [2];
        double** dx = new double* [punktyCalkowania * punktyCalkowania];
        double** dy = new double* [punktyCalkowania * punktyCalkowania];


        Jakobian[0] = new double[2];
        Jakobian[1] = new double[2];
        JakobianOdwrotny[0] = new double[2];
        JakobianOdwrotny[1] = new double[2];

        for (int i = 0; i < punktyCalkowania * punktyCalkowania; i++) {
            dx[i] = new double[4];
            dy[i] = new double[4];
        }

        for (int i = 0; i < nE; i++) {
            double w1 = 5.0 / 9.0, w2 = 8.0 / 9.0, w3 = w1, w4 = w1;
            for (int j = 0; j < punktyCalkowania * punktyCalkowania; j++) {
                //Macierz H
                obliczJakobian(i, j, pochodnaEta, pochodnaKsi, Jakobian);
                obliczJakobianOdwrotny(Jakobian, JakobianOdwrotny);
                obliczPochodneDxDy(dx, dy, pochodnaEta, pochodnaKsi, JakobianOdwrotny, j);
                obliczMacierzH(elements[i].MacierzH, dx, dy, Jakobian, j, w1 * w4);

                //Macierz C
                elements[i].obliczMacierzC(j, nc, specificHeat, density, Jakobian, punktyCalkowania, w1 * w4);
               
                
                //std::cout << "-----------DEBUG----------" << "\n";
                //std::cout << "Iteracja " << i << ", " << j << "\t";
                //std::cout << "w1 = " << w1 << "\t";
                //std::cout << "w2 = " << w2 << "\t";
                //std::cout << "w3 = " << w3 << "\t";
                //std::cout << "\n";
               
                w1 = w2;                
                w2 = w3;                
                w3 = w1;               

                if (j == 7) {
                    //środkowy element
                    w1 = 8.0 / 9.0;
                    w4 = w1;
                }                
            }

            mnożenieMacierzyH(elements[i].MacierzH, Jakobian, conductivity);
            //elements[i].printH();
        }
    }

    void obliczJakobian(int i, int j, double** dEta, double** dKsi, double** Jakobian) {
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                Jakobian[i][j] = 0;
            }
        }

        for (int k = 0; k < 4; k++) {
            Jakobian[0][0] += dKsi[j][k] * elements[i].nodes[k].y;
            //std::cout << Jakobian[0][0] << "\n";
            Jakobian[0][1] += dEta[j][k] * elements[i].nodes[k].y;
            //std::cout << Jakobian[0][1] << "\n";
            Jakobian[1][0] += dKsi[j][k] * elements[i].nodes[k].x;
            //std::cout << Jakobian[1][0] << "\n";
            Jakobian[1][1] += dEta[j][k] * elements[i].nodes[k].x;
            //std::cout << Jakobian[1][1] << "\n";
        }

        /*
        std::cout << "TEST JAKOBIAN:\n\n";
        std::cout << "[0][0] = "<<Jakobian[0][0] << "\n";
        std::cout << "[0][1] = "<<Jakobian[0][1] << "\n";
        std::cout << "[1][0] = "<<Jakobian[1][0] << "\n";
        std::cout << "[1][1] = "<<Jakobian[1][1] << "\n";
        std::cout << "\n\n";
        */
    }

    void obliczJakobianOdwrotny(double** Jakobian, double** JakobianOdwrotny) {
        double wyznacznik = (Jakobian[0][0] * Jakobian[1][1] - Jakobian[0][1] * Jakobian[1][0]);
        //std::cout << "Wyznacznik: " << wyznacznik << "\n";
        double wyznacznikOdwrotny = 1 / wyznacznik;
        //std::cout << "Wyznacznik odwrotny: " << wyznacznikOdwrotny << "\n";
        JakobianOdwrotny[0][0] = wyznacznikOdwrotny * Jakobian[0][0];
        JakobianOdwrotny[0][1] = wyznacznikOdwrotny * Jakobian[0][1];
        JakobianOdwrotny[1][0] = wyznacznikOdwrotny * Jakobian[1][0];
        JakobianOdwrotny[1][1] = wyznacznikOdwrotny * Jakobian[1][1];

        /*
        std::cout << "TEST JAKOBIAN ODWROTNY:\n\n";
        std::cout << JakobianOdwrotny[0][0] << "\n";
        std::cout << JakobianOdwrotny[0][1] << "\n";
        std::cout << JakobianOdwrotny[1][0] << "\n";
        std::cout << JakobianOdwrotny[1][1] << "\n";
        std::cout << "\n\n";
        */
    }

    void obliczPochodneDxDy(double** dx, double** dy, double** dEta, double** dKsi, double** JakobianOdwrotny, int i) {
        for (int j = 0; j < 4; j++) {
            dx[i][j] = JakobianOdwrotny[0][0] * dEta[i][j] + JakobianOdwrotny[0][1] * dKsi[i][j];
            dy[i][j] = JakobianOdwrotny[1][0] * dEta[i][j] + JakobianOdwrotny[1][1] * dKsi[i][j];
            //std::cout << "\n\nIteracja: (" << i << "," << j << ")\n dX = " << dX[i][j] << "\t" << "dY = " << dY[i][j] << "\n";
        }
    }

    void obliczMacierzH(double** h, double** dx, double** dy, double** J, int a, double w) {
        if (punktyCalkowania == 2) {
            //std::cout << "Dla pc = 2"\n;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    h[i][j] += (dx[a][i] * dx[a][j]);
                    h[i][j] += (dy[a][i] * dy[a][j]);
                    //std::cout << h[i][j] << "\n";
                }
            }
        }
        else if (punktyCalkowania == 3) {
            //std::cout << "Dla pc = 3"\n;
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    h[i][j] += (dx[a][i] * dx[a][j]) * w;
                    h[i][j] += (dy[a][i] * dy[a][j]) * w;
                    //std::cout << h[i][j] << "\n";
                }
            }
        }
    }

    void mnożenieMacierzyH(double** h, double** Jakobian, double t) {
        double wyznacznik = Jakobian[0][0] * Jakobian[1][1] - Jakobian[0][1] * Jakobian[1][0];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                h[i][j] *= wyznacznik;
                h[i][j] *= t;
                //std::cout << "t = " << t << "\n";
                //std::cout << h[i][j] << "\n";
            }
        }        
    }

    void agregacjaMacierzyH() {
        double** tabAgregacji = new double* [nN];
        for (int i = 0; i < nN; i++) {
            tabAgregacji[i] = new double[nN];
        }

        for (int i = 0; i < nE; i++) {            
            for (int i = 0; i < nN; i++) {                
                for (int j = 0; j < nN; j++) {
                    tabAgregacji[i][j] = 0;
                }                    
            }
            elements[i].agregacja(tabAgregacji);
            for (int i = 0; i < nN; i++) {
                for (int j = 0; j < nN; j++) {
                    globalneH[i][j] += tabAgregacji[i][j];
                }
            }
        }  


        //wyświetlGlobalneH();
    }

    void wyświetlGlobalneH() {
        //wyświetlenie po agregacji
        std::cout << "Globalna Macierz H: \n";
        for (int i = 0; i < nN; i++) {
            for (int j = 0; j < nN; j++) {
                std::cout << globalneH[i][j] << "\t";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    void wyświetlKońcoweH() {
        //wyświetlenie dla [H] = [H] + [C]/dT
        std::cout << "Koncowa Macierz H: \n";
        for (int i = 0; i < nN; i++) {
            for (int j = 0; j < nN; j++) {
                std::cout << globalneH[i][j] << "\t";
            }                
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    void dodajDoTablicyAgregacji(double** tab) {
        for (int i = 0; i < nN; i++) {
            for (int j = 0; j < nN; j++) {
                globalneH[i][j] += tab[i][j];
            }                
        }           
    }

    void dodajDoWektoraGlobalnego(double* wektor) {
        for (int i = 0; i < nN; i++) {
            globalnyP[i] += wektor[i];
        }            
    }

    void agregacjaMacierzyC() {
        double** tempTablicaAgregacji = new double* [nN];
        for (int i = 0; i < nN; i++) {
            tempTablicaAgregacji[i] = new double[nN];
        }

        for (int i = 0; i < nE; i++) {            
            for (int i = 0; i < nN; i++) {                
                for (int j = 0; j < nN; ++j) {
                    //Zerowanie wartości
                    tempTablicaAgregacji[i][j] = 0;
                }                    
            }
            elements[i].agregacjaC(tempTablicaAgregacji);
            for (int i = 0; i < nN; i++) {
                for (int j = 0; j < nN; j++) {
                    globalneC[i][j] += tempTablicaAgregacji[i][j];
                }
            }
        }
       //wyświetlGlobalneC();
    }

    void wyświetlGlobalneC() {
        std::cout << "Globalna Macierz C: " << "\n";
        for (int i = 0; i < nN; i++) {
            for (int j = 0; j < nN; j++) {
                std::cout << globalneC[i][j] << "\t";
            }                
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    void sumowanieMacierzyHorazC() {
        for (int i = 0; i < nN; i++) {
            for (int j = 0; j < nN; j++) {
                globalneH[i][j] += globalneC[i][j] / simulationStepTime;
            }                
        }            
        //wyświetlKońcoweH();
    }

    void wyświetlGlobalnyP() {
        std::cout << "Globalny wektor P: " << "\n";
        for (int i = 0; i < nN; i++) {
            std::cout << globalnyP[i] << "\t";
        }            
        std::cout << "\n\n";
    }

    void ostatnieObliczenieP() {
        //{P} + {[C]/dT * {T0}}
        //std::cout << "--------------------------------\n\n";
        double suma;        
        for (int i = 0; i < nN; i++) {
            suma = 0.0;
            for (int j = 0; j < nN; j++) {
                suma += (globalneC[i][j] / simulationStepTime) * T0[j];
            }
            koncowyP[i] = globalnyP[i] + suma;
            //std::cout << "Wektor koncowy:\t" << koncowyP[i] << "\n";
        }
        //std::cout << "--------------------------------\n\n";
    }
};

struct sciana {
    Element* w;
    double* n1;
    double* n2;
    double* n3;
    int rozmiar;

    sciana() {};
    //Uzupełnienie ściany dla pc == 2
    sciana(int size, double ksi1, double eta1, double ksi2, double eta2) {
        rozmiar = size;
        w = new Element[size];
        n1 = new double[4];
        n2 = new double[4];
        n3 = new double[4];

        n1[0] = 0.25 * (1.0 - ksi1) * (1 - eta1);
        n1[1] = 0.25 * (1 + ksi1) * (1 - eta1);
        n1[2] = 0.25 * (1 + ksi1) * (1 + eta1);
        n1[3] = 0.25 * (1 - ksi1) * (1 + eta1);

        n2[0] = 0.25 * (1.0 - ksi2) * (1 - eta2);
        n2[1] = 0.25 * (1 + ksi2) * (1 - eta2);
        n2[2] = 0.25 * (1 + ksi2) * (1 + eta2);
        n2[3] = 0.25 * (1 - ksi2) * (1 + eta2);

        n3[0] = 0;
        n3[1] = 0;
        n3[2] = 0;
        n3[3] = 0;
    }

    //Uzupełnienie ściany dla pc == 3
    sciana(int size, double ksi1, double eta1, double ksi2, double eta2, double ksi3, double eta3) {
        rozmiar = size;
        w = new Element[size];
        n1 = new double[4];
        n2 = new double[4];
        n3 = new double[4];

        n1[0] = 0.25 * (1.0 - ksi1) * (1 - eta1);
        n1[1] = 0.25 * (1 + ksi1) * (1 - eta1);
        n1[2] = 0.25 * (1 + ksi1) * (1 + eta1);
        n1[3] = 0.25 * (1 - ksi1) * (1 + eta1);

        n2[0] = 0.25 * (1.0 - ksi2) * (1 - eta2);
        n2[1] = 0.25 * (1 + ksi2) * (1 - eta2);
        n2[2] = 0.25 * (1 + ksi2) * (1 + eta2);
        n2[3] = 0.25 * (1 - ksi2) * (1 + eta2);

        n3[0] = 0.25 * (1.0 - ksi3) * (1 - eta3);
        n3[1] = 0.25 * (1 + ksi3) * (1 - eta3);
        n3[2] = 0.25 * (1 + ksi3) * (1 + eta3);
        n3[3] = 0.25 * (1 - ksi3) * (1 + eta3);
    }
};

struct Element2D_4 {
    double punktCalkowania = 1.0 / sqrt(3.0);
    double** dKsi;
    double** dEta;
    double** nc;

    sciana lewaSciana, prawaSciana, gornaSciana, dolnaSciana;
    Grid* siatka;

    Element2D_4(Grid* siatka) {
        this->siatka = siatka;

        int pc = siatka->punktyCalkowania; //zmienna pomocnicza
        //Podział w zależności od ilości punktów całkowania

        dKsi = new double* [pc * pc];
        dEta = new double* [pc * pc];
        nc = new double* [pc * pc];

        for (int i = 0; i < pc * pc; i++) {
            //Końcowe wymiary - 4x4 dla pc == 2 lub 9x4 dla pc ==3

            dKsi[i] = new double[4];
            dEta[i] = new double[4];
            nc[i] = new double[4];
        }

        //Uzupełnianie wartości do tablic
        if (pc == 2) {
            uzupelnijFunkcjeKsztaltu(nc,-(1.0/sqrt(3.0)),-(1.0/sqrt(3.0)));
        }
        else if (pc == 3) {
            uzupelnijFunkcjeKsztaltu(nc, sqrt(3.0 / 5.0), 1.0);
        }
        siatka->nc = nc;
        
        //Przypisanie punktów do ścian
        if (pc == 2) {
            lewaSciana = sciana(siatka->nH, -1.0, (1.0/sqrt(3.0)), -1.0, -(1.0/sqrt(3.0)));
            prawaSciana = sciana(siatka->nH, 1.0, -(1.0/sqrt(3.0)), 1.0, (1.0/sqrt(3.0)));
            gornaSciana = sciana(siatka->nB, -(1.0/sqrt(3.0)), 1.0, (1.0/sqrt(3.0)), 1.0);
            dolnaSciana = sciana(siatka->nB, -(1.0/sqrt(3.0)), -1.0, (1.0/sqrt(3.0)), -1.0);
        }

        else if (pc == 3) {
            lewaSciana = sciana(siatka->nH, -1.0, (sqrt(3.0/5.0)), -1.0, 0, -1.0, -(sqrt(3.0/5.0)));
            prawaSciana = sciana(siatka->nH, 1.0, -(sqrt(3.0/5.0)), 1.0, 0, 1.0, (sqrt(3.0/5.0)));
            gornaSciana = sciana(siatka->nB, (sqrt(3.0/5.0)), 1.0, 0, 1.0, -(sqrt(3.0/5.0)), 1.0);
            dolnaSciana = sciana(siatka->nB, -(sqrt(3.0/5.0)), -1.0, 0, -1.0, (sqrt(3.0/5.0)), -1.0);
        }
        uzupelnijSciany(siatka->elements);
    }

    void uzupelnijFunkcjeKsztaltu(double** tablica, double ksi, double eta) {
        //Podział w zależności od ilości punktów całkowania

        if (siatka->punktyCalkowania == 2) {
            tablica[0][0] = 0.25 * (1.0 - ksi) * (1 - eta);
            tablica[0][1] = 0.25 * (1 + ksi) * (1 - eta);
            tablica[0][2] = 0.25 * (1 + ksi) * (1 + eta);
            tablica[0][3] = 0.25 * (1 - ksi) * (1 + eta);

            ksi = -ksi;
            tablica[1][0] = 0.25 * (1.0 - ksi) * (1 - eta);
            tablica[1][1] = 0.25 * (1 + ksi) * (1 - eta);
            tablica[1][2] = 0.25 * (1 + ksi) * (1 + eta);
            tablica[1][3] = 0.25 * (1 - ksi) * (1 + eta);

            eta = -eta;
            tablica[2][0] = 0.25 * (1.0 - ksi) * (1 - eta);
            tablica[2][1] = 0.25 * (1 + ksi) * (1 - eta);
            tablica[2][2] = 0.25 * (1 + ksi) * (1 + eta);
            tablica[2][3] = 0.25 * (1 - ksi) * (1 + eta);

            ksi = -ksi;
            tablica[3][0] = 0.25 * (1.0 - ksi) * (1 - eta);
            tablica[3][1] = 0.25 * (1 + ksi) * (1 - eta);
            tablica[3][2] = 0.25 * (1 + ksi) * (1 + eta);
            tablica[3][3] = 0.25 * (1 - ksi) * (1 + eta);
        }
        else if (siatka->punktyCalkowania == 3) {
            //zmienne pomocnicze
            double a = ksi;

            ksi = -a;
            eta = -a;
            tablica[0][0] = 0.25 * (1.0 - ksi) * (1 - eta);
            tablica[0][1] = 0.25 * (1 + ksi) * (1 - eta);
            tablica[0][2] = 0.25 * (1 + ksi) * (1 + eta);
            tablica[0][3] = 0.25 * (1 - ksi) * (1 + eta);

            ksi = 0;
            tablica[1][0] = 0.25 * (1.0 - ksi) * (1 - eta);
            tablica[1][1] = 0.25 * (1 + ksi) * (1 - eta);
            tablica[1][2] = 0.25 * (1 + ksi) * (1 + eta);
            tablica[1][3] = 0.25 * (1 - ksi) * (1 + eta);

            ksi = a;
            tablica[2][0] = 0.25 * (1.0 - ksi) * (1 - eta);
            tablica[2][1] = 0.25 * (1 + ksi) * (1 - eta);
            tablica[2][2] = 0.25 * (1 + ksi) * (1 + eta);
            tablica[2][3] = 0.25 * (1 - ksi) * (1 + eta);

            eta = 0;
            tablica[3][0] = 0.25 * (1.0 - ksi) * (1 - eta);
            tablica[3][1] = 0.25 * (1 + ksi) * (1 - eta);
            tablica[3][2] = 0.25 * (1 + ksi) * (1 + eta);
            tablica[3][3] = 0.25 * (1 - ksi) * (1 + eta);

            eta = a;
            tablica[4][0] = 0.25 * (1.0 - ksi) * (1 - eta);
            tablica[4][1] = 0.25 * (1 + ksi) * (1 - eta);
            tablica[4][2] = 0.25 * (1 + ksi) * (1 + eta);
            tablica[4][3] = 0.25 * (1 - ksi) * (1 + eta);

            ksi = 0;
            tablica[5][0] = 0.25 * (1.0 - ksi) * (1 - eta);
            tablica[5][1] = 0.25 * (1 + ksi) * (1 - eta);
            tablica[5][2] = 0.25 * (1 + ksi) * (1 + eta);
            tablica[5][3] = 0.25 * (1 - ksi) * (1 + eta);

            ksi = -a;
            tablica[6][0] = 0.25 * (1.0 - ksi) * (1 - eta);
            tablica[6][1] = 0.25 * (1 + ksi) * (1 - eta);
            tablica[6][2] = 0.25 * (1 + ksi) * (1 + eta);
            tablica[6][3] = 0.25 * (1 - ksi) * (1 + eta);

            eta = 0;
            tablica[7][0] = 0.25 * (1.0 - ksi) * (1 - eta);
            tablica[7][1] = 0.25 * (1 + ksi) * (1 - eta);
            tablica[7][2] = 0.25 * (1 + ksi) * (1 + eta);
            tablica[7][3] = 0.25 * (1 - ksi) * (1 + eta);

            ksi = 0;
            tablica[8][0] = 0.25 * (1.0 - ksi) * (1 - eta);
            tablica[8][1] = 0.25 * (1 + ksi) * (1 - eta);
            tablica[8][2] = 0.25 * (1 + ksi) * (1 + eta);
            tablica[8][3] = 0.25 * (1 - ksi) * (1 + eta);
        }

    }

    void uzupelnijSciany(Element* el) {
        //zmienne pomocnicze
        int lewa = 0, prawa = 0, gora = 0, dol = 0;
        for (int i = 0; i < siatka->nE; i++) {
            if ((el[i].nodes[1].warunekBrzegowy == 1) && (el[i].nodes[2].warunekBrzegowy == 1)) {
                prawaSciana.w[prawa] = el[i];
                prawa++;
            }
            if ((el[i].nodes[3].warunekBrzegowy == 1) && (el[i].nodes[0].warunekBrzegowy == 1)) {
                lewaSciana.w[lewa] = el[i];
                lewa++;
            }
            if ((el[i].nodes[0].warunekBrzegowy == 1) && (el[i].nodes[1].warunekBrzegowy == 1)) {
                dolnaSciana.w[dol] = el[i];
                dol++;
            }
            if ((el[i].nodes[2].warunekBrzegowy == 1) && (el[i].nodes[3].warunekBrzegowy == 1)) {
                gornaSciana.w[gora] = el[i];
                gora++;
            }
        }
    }

    void wyswietlPochodneKsiorazEta() {

        if (siatka->punktyCalkowania == 2) {

            std::cout << "\ndN/dKsi:\n";
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    std::cout << dKsi[i][j] << "\t";
                }
                std::cout << "\n";
            }

            std::cout << "\ndN/dEta:\n";
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    std::cout << dEta[i][j] << "\t";
                }                    
                std::cout << "\n";
            }
        }

        else if (siatka->punktyCalkowania == 3) {

            std::cout << "\ndN/dKsi:\n";
            for (int i = 0; i < 9; i++) {
                for (int j = 0; j < 4; j++) {
                    std::cout << dKsi[i][j] << "\t";
                }                    
                std::cout << "\n";
            }

            std::cout << "\ndN/dEta:\n";
            for (int i = 0; i < 9; i++) {
                for (int j = 0; j < 4; j++) {
                    std::cout << dEta[i][j] << "\t";
                }                    
                std::cout << "\n";
            }
        }
    }


    double fun1(double x) {
        return (1.0 / 4.0) * (1 - x);
    }

    double fun2(double x) {
        return (1.0 / 4.0) * (1 + x);
    }

    double fun3(double x) {
        return (1.0 / 4.0) * (x - 1);
    }

    void obliczPochodneKsiEta(int pc) {
        if (pc == 2){
            double x = 1.0 / sqrt(3.0);

            double tabXY[4][2] = { {-x,-x},
                                    {x,-x},
                                    {x,x}, 
                                    {-x, x} };


            for (int i = 0; i < 4; i++) {
                dKsi[i][0] = fun3(tabXY[i][0]);
                dKsi[i][1] = -fun2(tabXY[i][0]);
                dKsi[i][2] = fun2(tabXY[i][0]);
                dKsi[i][3] = fun1(tabXY[i][0]);

                dEta[i][0] = fun3(tabXY[i][1]);
                dEta[i][1] = fun1(tabXY[i][1]);
                dEta[i][2] = fun2(tabXY[i][1]);
                dEta[i][3] = -fun2(tabXY[i][1]);
            }
        }
        else if (pc == 3) {
            double x = sqrt(3.0 / 5.0), x1 = 0.0;
            double tabXY[9][2] = { {-x,-x},
                                    {x1,-x},
                                    {x,-x},
                                    {x, x1},
                                    {x,x},
                                    {x1,x},
                                    {-x,x},
                                    {-x,x1},
                                    {x1,x1} };


            for (int i = 0; i < 9; i++) {
                dKsi[i][0] = fun3(tabXY[i][0]);
                dKsi[i][1] = -fun2(tabXY[i][0]);
                dKsi[i][2] = fun2(tabXY[i][0]);
                dKsi[i][3] = fun1(tabXY[i][0]);

                dEta[i][0] = fun3(tabXY[i][1]);
                dEta[i][1] = fun1(tabXY[i][1]);
                dEta[i][2] = fun2(tabXY[i][1]);
                dEta[i][3] = -fun2(tabXY[i][1]);
            }
        }
    }

    void obliczHbc() {
        for (int i = 0; i < siatka->nB - 1; i++) {
            uzupełnijHbc(dolnaSciana.w[i], dolnaSciana, 0);
            uzupełnijHbc(gornaSciana.w[i], gornaSciana, 2);
        }
        for (int i = 0; i < siatka->nH - 1; i++) {
            uzupełnijHbc(lewaSciana.w[i], lewaSciana, 3);
            uzupełnijHbc(prawaSciana.w[i], prawaSciana, 1);
        }

        /*
        for(int i=0;i<siatka->nE;++i)
        {
            std::cout<<"HBC"<< "\n";
            siatka->elements[i].wyświetlHBC();
        }
        //Wyświetlanie HBC


        for(int i=0; i<siatka->nB-1; ++i)
        {
            std::cout<<"HBC"<< "\n";
            dolnaSciana.w[i].wyświetlHBC();
            std::cout<< "\n";
            std::cout<<"HBC"<< "\n";
            gornaSciana.w[i].wyświetlHBC();
            std::cout<< "\n";
        }

        for(int i=0; i<siatka->nH-1; ++i)
        {
            std::cout<<"HBC"<< "\n";
            lewaSciana.w[i].wyświetlHBC();
            std::cout<< "\n";
            std::cout<<"HBC"<< "\n";
            prawaSciana.w[i].wyświetlHBC();
            std::cout<< "\n";
        }
        */
    }

    void uzupełnijHbc(Element e, sciana wall, int surfaceId) {
        double det = e.długość(surfaceId);
        //e.MacierzC[0][0] = 0;
        //e.MacierzHbc[0][0] = 0;
        //double det = 0.00015625;
        if (siatka->punktyCalkowania == 2) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    e.MacierzHbc[j][k] += siatka->alfa * det / 2 * wall.n1[j] * wall.n1[k];
                    //std::cout << "Test dla mianownika: (" << j <<", " << k << ") = " << (2 * wall.n1[j] * wall.n1[k]);
                    //std::cout << e.MacierzHbc[j][k] << "\n";
                    e.MacierzHbc[j][k] += siatka->alfa * det / 2 * wall.n2[j] * wall.n2[k];
                    //std::cout << e.MacierzHbc[j][k] << "\n";
                }
                //std::cout<<"\n";
            }                
        }
        else if (siatka->punktyCalkowania == 3) {
            double w1 = 5.0 / 9.0, w2 = 8.0 / 9.0;
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    e.MacierzHbc[j][k] += siatka->alfa * det / 2 * w1 * wall.n1[j] * wall.n1[k];
                    e.MacierzHbc[j][k] += siatka->alfa * det / 2 * w2 * wall.n2[j] * wall.n2[k];
                    e.MacierzHbc[j][k] += siatka->alfa * det / 2 * w1 * wall.n3[j] * wall.n3[k];
                }
            }                
        }
    }



    void agregacjaWektoraP() {

        obliczWektoryP();
        double* wektorTemp = new double[siatka->nN];
        for (int i = 0; i < siatka->nE; i++) {            
            for (int i = 0; i < siatka->nN; i++) {
                wektorTemp[i] = 0;
            }
            siatka->elements[i].agregacjaP(wektorTemp);
            siatka->dodajDoWektoraGlobalnego(wektorTemp);
        }
        //siatka->wyświetlGlobalnyP();
    }

    

    void obliczWektoryP() {
        for (int i = 0; i < siatka->nB - 1; i++) {
            dolnaSciana.w[i].obliczanieP(dolnaSciana.n1, dolnaSciana.n2, dolnaSciana.n3, siatka->ambientTemperature, 0, siatka->punktyCalkowania);
            gornaSciana.w[i].obliczanieP(gornaSciana.n1, gornaSciana.n2, gornaSciana.n3, siatka->ambientTemperature, 2, siatka->punktyCalkowania);
        }
        for (int i = 0; i < siatka->nH - 1; i++) {
            lewaSciana.w[i].obliczanieP(lewaSciana.n1, lewaSciana.n2, lewaSciana.n3, siatka->ambientTemperature, 3, siatka->punktyCalkowania);
            prawaSciana.w[i].obliczanieP(prawaSciana.n1, prawaSciana.n2, prawaSciana.n3, siatka->ambientTemperature, 1, siatka->punktyCalkowania);
        }
    }

    void wynik() {
        double* rozw = rozwiążRównanie(siatka->globalneH, siatka->koncowyP, siatka->nN);
        //std::cout << "Test case : \t" << *sol << "\n";
        double min = rozw[0], max = 0.0;
        for (int i = 0; i < siatka->nN; i++) {
            siatka->T0[i] = rozw[i];
            siatka->koncowyP[i] = 0;
            //std::cout << sol[i] << "\n";
            if (rozw[i] < min) {
                min = rozw[i];
            }                
            else if (rozw[i] > max) {
                max = rozw[i];
            }                
        }
        std::cout << "MinTemp: " << min << " MaxTemp: " << max << "\n";
        std::cout << "\n";

    }

    double* rozwiążRównanie(double** a, double* b, int num) {
        double* n = new double[num];
        double* x1 = new double[num];
        double* x2 = new double[num];
        double** M = new double* [num];

        for (int i = 0; i < num; i++) {
            x1[i] = 0.0;
            x2[i] = 0.0;
            n[i] = 0.0;
        }

        for (int i = 0; i < num; i++) {
            M[i] = new double[num];
            for (int j = 0; j < num; ++j) {
                M[i][j] = 0.0;
            }                
        }

        for (int i = 0; i < num; i++) {
            n[i] = 1 / a[i][i];
        }
        // Calculate M = -D^-1 (L + U)
        for (int i = 0; i < num; i++) {
            for (int j = 0; j < num; j++) {
                if (i == j) {
                    M[i][j] = 0.0;
                }
                else {
                    M[i][j] = -(a[i][j] * n[i]);
                }                    
            }
        }           

        for (int k = 0; k < 100; k++) {
            for (int i = 0; i < num; i++) {
                x2[i] = n[i] * b[i];
                for (int j = 0; j < num; j++) {
                    x2[i] += M[i][j] * x1[j];
                }                    
            }
            for (int i = 0; i < num; i++) {
                x1[i] = x2[i];
            }                
        }
        return x1;
    }
};

int main()
{
    std::cout << "Obliczanie temperatur dla procesu niestacjonarnego przepływu ciepła\n";

    std::cout.precision(3);
    std::cout<<std::fixed;


    /*
    double h = 0.2;
    double b = 0.1;
    int nH = 7;
    int nB = 5;
    int nN = nH * nB;
    int nE = (nH - 1) * (nB - 1);    
    */    
    //---------------------------------------------------------
    //Zmienić wartości dla innych przykładów
    //Dla pierwszego przypadku - siatka 4x4
    
    int tempPocz = 100;
    int czasSymulacji = 500;
    int krokCzasowy = 50;
    int tempOtocz = 1200;
    int alfa = 300;
    int specificHeat = 700;
    int conductivity = 25;
    int density = 7800;

    double h = 0.1, b = 0.1;
    int nH = 4, nB = 4;
    int punktyCalkowania = 2; //Zmienić na 3 dla trójpunktowego
    
    //---------------------------------------------------------

    //Dla drugiego przypadku - siatka 31x31
    /*
    int tempPocz = 100;
    int czasSymulacji = 100; //Można zmienić na jakiś mniejszy, by program szybciej kończył działanie
    int krokCzasowy = 1;
    int tempOtocz = 1200;
    int alfa = 300;
    int specificHeat = 700;
    int conductivity = 25;
    int density = 7800;

    double h = 0.1, b = 0.1;
    int nH = 31, nB = 31;
    int punktyCalkowania = 2; //Zmienić na 3 dla trójpunktowego
    */


    Grid siatka(h, b, nH, nB,punktyCalkowania,tempPocz,krokCzasowy,tempOtocz,alfa,specificHeat,conductivity,density);
    siatka.wyswietlSiatke();
    // 
    // 
    //siatka.show();
    
    Element2D_4 elementWzorcowy(&siatka);

    
    elementWzorcowy.obliczPochodneKsiEta(punktyCalkowania);
    //elementWzorcowy.wyswietlPochodneKsiorazEta();


    siatka.funkcjaJakobian(elementWzorcowy.dKsi, elementWzorcowy.dEta);
    
    elementWzorcowy.obliczHbc();
    
    siatka.agregacjaMacierzyH();
   
    elementWzorcowy.agregacjaWektoraP();

    siatka.agregacjaMacierzyC();

    siatka.sumowanieMacierzyHorazC();

    for (int i = 0; i < czasSymulacji / krokCzasowy; i++) {
        std::cout << "Iteracja "<< i+1 << "\n";
        siatka.ostatnieObliczenieP();
        elementWzorcowy.wynik();
    }
    
}

// Uruchomienie programu: Ctrl + F5 lub menu Debugowanie > Uruchom bez debugowania
// Debugowanie programu: F5 lub menu Debugowanie > Rozpocznij debugowanie

// Porady dotyczące rozpoczynania pracy:
//   1. Użyj okna Eksploratora rozwiązań, aby dodać pliki i zarządzać nimi
//   2. Użyj okna programu Team Explorer, aby nawiązać połączenie z kontrolą źródła
//   3. Użyj okna Dane wyjściowe, aby sprawdzić dane wyjściowe kompilacji i inne komunikaty
//   4. Użyj okna Lista błędów, aby zobaczyć błędy
//   5. Wybierz pozycję Projekt > Dodaj nowy element, aby utworzyć nowe pliki kodu, lub wybierz pozycję Projekt > Dodaj istniejący element, aby dodać istniejące pliku kodu do projektu
//   6. Aby w przyszłości ponownie otworzyć ten projekt, przejdź do pozycji Plik > Otwórz > Projekt i wybierz plik sln

/*

    void obliczMacierzH(double** H, double** dX, double** dY, double** Jakobian) {

        //sprawdzic, gdzie tutaj konkretnie jest jakis error (przy mnozeniu macierzy razy ta sama macierz transponowana
        //Porownać wyniki w programie z wynikami z kartki
        //Obliczamy wartosci w nawiasach dla wzoru z macierza H


        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                for (int k = 0; k < 4; k++) {
                    H[j][k] += dX[i][j] * dX[i][k];
                    //std::cout << H[j][k] << "\t";

                    H[j][k] += dY[i][j] * dY[i][k];
                }
                //std::cout << "\n";
            }
            //std::cout << "\n";
        }

        //test

        /*std::cout << "----------\nWYSWIETLAM MACIERZ H\n------------\n";
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                std::cout << H[i][j] << "\t";
            }
            std::cout << "\n";
        }

        //Do wartosci w nawiasach dla kazdego punktu mnozymy jeszcze k(t) i dV
double dV = (Jakobian[0][0] * Jakobian[1][1]) - (Jakobian[0][1] * Jakobian[1][0]);
std::cout << "Wyznacznik i dV - " << dV << "\n";
for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
        H[i][j] *= 25.0;
        H[i][j] *= dV;
    }
}
    }
    
    */