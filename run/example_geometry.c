void example_geometry() {
    double mm = 0.1;
    double degree = 1;
    cout << "Root, TShield example" << endl;
    gSystem->Load("libGeom.so");
    new TGeoManager("World", "Geometry for example.c root file.");

    TGeoElementTable *table = gGeoManager->GetElementTable();
    TGeoMaterial *he = new TGeoMaterial("He", table->GetElement(2), 1);
    TGeoMaterial *al = new TGeoMaterial("Al", table->GetElement(13), 1);
    TGeoMaterial *cl = new TGeoMaterial("Cl", table->GetElement(17), 1);
    TGeoMaterial *fe = new TGeoMaterial("Fe", table->GetElement(26), 1);
    TGeoMaterial *h  = new TGeoMaterial("H" , table->GetElement(1), 1);
    TGeoMaterial *pb = new TGeoMaterial("Pb", table->GetElement(82), 1);
    TGeoMaterial *vac = new TGeoMaterial("Vacuum", 0, 0, 0);

    TGeoMedium *mvac = new TGeoMedium("Vacuum", 1, vac);
    TGeoMedium *mal = new TGeoMedium("Al", 1, al);
    TGeoMedium *mhe = new TGeoMedium("He", 1, he);
    TGeoMedium *mcl = new TGeoMedium("Ar", 1, cl);
    TGeoMedium *mfe = new TGeoMedium("Fe", 1, fe);
    TGeoMedium *mh  = new TGeoMedium("H" , 1, h);
    TGeoMedium *mpb = new TGeoMedium("Pb", 1, pb);
    TGeoMixture *mhehar = new TGeoMixture("He+H+Ar", 1, 0.5);
    mhehar->AddElement(he, 0.5);
    mhehar->AddElement(h, 0.1);
    mhehar->AddElement(cl, 0.4);

    TGeoVolume *top = gGeoManager->MakeBox("Top", mfe, 100 * mm, 50 * mm, 50 * mm);

    TGeoPgon *pcon = new TGeoPcon(0, 360, 5);
    pcon->DefineSection(0, 0, 0, 40 * mm);
    pcon->DefineSection(1, 20 * mm, 0, 40 * mm);
    pcon->DefineSection(2, 30 * mm, 0, 30 * mm);
    pcon->DefineSection(3, 40 * mm, 0, 40 * mm);
    pcon->DefineSection(4, 100 * mm, 0, 40 * mm);
    TGeoVolume *vpcon = new TGeoVolume("PCon", pcon, mal);
    TGeoSphere *sph = new TGeoSphere(0, 30 * mm, 0, 180 * degree, 0, 360 * degree);
    TGeoVolume *vsph = new TGeoVolume("Sph", sph, mcl);
    TGeoTube *tub = new TGeoTube(5 * mm, 6 * mm, 20 * mm);
    TGeoVolume *vtub = new TGeoVolume("Tub", tub, mpb);
    TGeoTube *tubIn = new TGeoTube(0 * mm, 5 * mm, 20 * mm);
    TGeoVolume *vtubIn = new TGeoVolume("TubIn", tubIn, mcl);
    TGeoConeSeg *con = new TGeoConeSeg(10 * mm, 0, 5 * mm, 0, 40 * mm , 45 * degree, (360 - 45)*degree);
//     TGeoCone* con = new TGeoCone(10*mm,0, 5 * mm, 0, 40*mm);
    TGeoVolume *vcon = new TGeoVolume("Con", con, mal);
    TGeoBBox *box = new TGeoBBox(10 * mm, 10 * mm, 10 * mm);
    TGeoVolume *vbox = new TGeoVolume("Box", box, mal);
    TGeoPgon *pgon = new TGeoPgon(0, 360, 3, 2);
    pgon->DefineSection(0, 0, 0, 10 * mm);
    pgon->DefineSection(1, 40 * mm, 0, 25 * mm);
    TGeoVolume *vpgon = new TGeoVolume("PGon", pgon, mal);

    TGeoHMatrix *mpcon1 = new TGeoHMatrix(); mpcon1->SetDz(20*mm);
    TGeoHMatrix *mpcon2 = new TGeoHMatrix(); mpcon2->SetDz(70*mm);
    vtub->AddNode(vtubIn,0);
    vpcon->AddNode(vtub,0,mpcon1);
    vpcon->AddNode(vsph,1,mpcon2);
    TGeoHMatrix *mtop1 = new TGeoHMatrix(); mtop1->RotateY(-90); mtop1->SetDx(100* mm);
    TGeoHMatrix *mtop2 = new TGeoHMatrix(); mtop2->RotateZ(180); mtop2->RotateY(90); mtop2->SetDx(-10 * mm);
    TGeoHMatrix *mtop3 = new TGeoHMatrix(); mtop3->SetDx(-30* mm);
    TGeoHMatrix *mtop4 = new TGeoHMatrix(); mtop4->SetDx(-50* mm);
    TGeoHMatrix *mtop5 = new TGeoHMatrix(); mtop5->RotateY(-90); mtop5->SetDx(-60* mm);
    top->AddNode(vpcon,0,mtop1);
    top->AddNode(vcon,1,mtop2);
    top->AddNode(vbox,2,mtop3);
    top->AddNode(vbox,2,mtop4);
   top->AddNode(vpgon,0,mtop5);
    

    gGeoManager->SetTopVolume(top);
    gGeoManager->CloseGeometry();
    top->SetLineColor(kMagenta);
    gGeoManager->SetTopVisible();

    gGeoManager->GetTopVolume()->Draw();

    gSystem->Load("libTShield.so");
//    TShield b;
//    b.SetGeometry(gGeoManager);
//    b.PrintGeometry();
}
