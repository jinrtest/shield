void printElement(Element element) {
    printf("Element:\t{\t%f,\t%f,\t%f,\t%f,\t%f,\t%f,\t%f}\n",
           element.Nuclid, element.Conc, element.Density,
           element.Z, element.A, element.PureDensity, element.ionEv);
}
void printMediumData(MediumData medium) {
    printf("Medium:\n");
    printf("\tType: %i\n", medium.nType);
    printf("\tRho: %f\n", medium.Rho);
    printf("\tCount of elements: %i\n", medium.nChemEl);
    for (int i = 0; i < 24; i++) {
        printf("\tElement %i\n\t", i);
        printElement(medium.Elements[i]);
    }
    printf("\n");
}
void printSGeoBody(SGeoBody body) {
    printf("Body: %i\n", body.type);
    printf("\t");
    for (int i = 0; i < 36; i++) {
        printf("%f ", body.parameters[i]);
    }
    printf("\n");
}
void printSGeoZone(SGeoZone zone) {
    printf("Zone:\n");
    printf("\tNumber of medium: %i\n", zone.mediaNumber);
    printf("\tCount of elements: %i\n", zone.countELements);
    printf("\t");
    for (int i = 0; i < zone.countELements; i++) {
        printf("%i ", zone.definition[2 * i]*zone.definition[2 * i + 1]);
    }
    printf("\n");
}
