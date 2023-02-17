/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "faceProbes.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "IOmanip.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceProbes, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        faceProbes,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::faceProbes::findElements(const fvMesh& mesh)
{
    if (debug)
    {
        Info<< "faceProbes: resetting sample locations" << endl;
    }

    elementList_.clear();
    elementList_.setSize(size());

    faceList_.clear();
    faceList_.setSize(size());

    forAll(*this, faceProbei)
    {
        const vector& location = operator[](faceProbei);
	vector& fluxDir = fluxVector_;//{1, 0, 0};
        const label celli = mesh.findCell(location);
	labelList fluxFaces;
	fluxFaces.setSize(2);
        elementList_[faceProbei] = celli;

        if (celli != -1)
        {
            const labelList& cellFaces = mesh.cells()[celli];
	    const vectorField& faceCentres = mesh.faceCentres();
	    const surfaceVectorField& surfaces = mesh.Sf();
            label choosenFaceID = cellFaces[0];
	    double minDotProduct1 = 1.0, minDotProduct2 = 1.0;

	    forAll(cellFaces, i)
            {
                label facei = cellFaces[i];
                vector surface = surfaces[facei]/mag(surfaces[facei]);
		double dotProduct = fluxDir & surface;
		double diffDotProduct = mag(mag(dotProduct) - 1.0); 
		//Debug
		/*Info << "Face ID: " << facei << " Surface vector: " << surface <<
			" Dot Product: " << diffDotProduct << " face centre: " <<
		       	faceCentres[facei] << endl;
	        */
		if (diffDotProduct < minDotProduct1)
        	{
                    minDotProduct1 = diffDotProduct;
                    fluxFaces[0] = facei;
                }

		else if (diffDotProduct < minDotProduct2)
		{
                    minDotProduct2 = diffDotProduct;
                    fluxFaces[1] = facei;
		}

	    }

	    if (
	        (fluxDir & faceCentres[fluxFaces[0]]) < (fluxDir & faceCentres[fluxFaces[1]])
		and mesh.isInternalFace(fluxFaces[0])		
   	       )
	    {
	        choosenFaceID = fluxFaces[0];
	    }
	    else if(mesh.isInternalFace(fluxFaces[1]))
	    {
	        choosenFaceID = fluxFaces[1];
	    }
	    else 
	    {
		choosenFaceID = fluxFaces[0];
	    }
	    
	    //Info << "Choosen ID: " << choosenFaceID << endl;
            faceList_[faceProbei] = choosenFaceID;
        }
        else
        {
            faceList_[faceProbei] = -1;
        }

        if (debug && (elementList_[faceProbei] != -1 || faceList_[faceProbei] != -1))
        {
            Pout<< "faceProbes : found point " << location
                << " in cell " << elementList_[faceProbei]
                << " and face " << faceList_[faceProbei] << endl;
        }
    }


    // Check if all faceProbes have been found.
    forAll(elementList_, faceProbei)
    {
        const vector& location = operator[](faceProbei);
        label celli = elementList_[faceProbei];
        label facei = faceList_[faceProbei];

        // Check at least one processor with cell.
        reduce(celli, maxOp<label>());
        reduce(facei, maxOp<label>());

        if (celli == -1)
        {
            if (Pstream::master())
            {
                WarningInFunction
                    << "Did not find location " << location
                    << " in any cell. Skipping location." << endl;
            }
        }
        else if (facei == -1)
        {
            if (Pstream::master())
            {
                WarningInFunction
                    << "Did not find location " << location
                    << " in any face. Skipping location." << endl;
            }
        }
        else
        {
            // Make sure location not on two domains.
            if (elementList_[faceProbei] != -1 && elementList_[faceProbei] != celli)
            {
                WarningInFunction
                    << "Location " << location
                    << " seems to be on multiple domains:"
                    << " cell " << elementList_[faceProbei]
                    << " on my domain " << Pstream::myProcNo()
                        << " and cell " << celli << " on some other domain."
                    << endl
                    << "This might happen if the faceProbe location is on"
                    << " a processor patch. Change the location slightly"
                    << " to prevent this." << endl;
            }

            if (faceList_[faceProbei] != -1 && faceList_[faceProbei] != facei)
            {
                WarningInFunction
                    << "Location " << location
                    << " seems to be on multiple domains:"
                    << " cell " << faceList_[faceProbei]
                    << " on my domain " << Pstream::myProcNo()
                        << " and face " << facei << " on some other domain."
                    << endl
                    << "This might happen if the faceProbe location is on"
                    << " a processor patch. Change the location slightly"
                    << " to prevent this." << endl;
            }
        }
    }
}

Foam::label Foam::faceProbes::prepare()
{
    const label nFields = classifyFields();

    // adjust file streams
    if (Pstream::master())
    {
        wordHashSet currentFields;

        currentFields.insert(scalarFields_);
        currentFields.insert(vectorFields_);
        currentFields.insert(sphericalTensorFields_);
        currentFields.insert(symmTensorFields_);
        currentFields.insert(tensorFields_);

        currentFields.insert(surfaceScalarFields_);
        currentFields.insert(surfaceVectorFields_);
        currentFields.insert(surfaceSphericalTensorFields_);
        currentFields.insert(surfaceSymmTensorFields_);
        currentFields.insert(surfaceTensorFields_);

        if (debug)
        {
            Info<< "Probing fields: " << currentFields << nl
                << "Probing locations: " << *this << nl
                << endl;
        }


        fileName faceProbeDir;
        fileName faceProbeSubDir = name();

        if (mesh_.name() != polyMesh::defaultRegion)
        {
            faceProbeSubDir = faceProbeSubDir/mesh_.name();
        }
        faceProbeSubDir = "postProcessing"/faceProbeSubDir/mesh_.time().timeName();

        if (Pstream::parRun())
        {
            // Put in undecomposed case
            // (Note: gives problems for distributed data running)
            faceProbeDir = mesh_.time().path()/".."/faceProbeSubDir;
        }
        else
        {
            faceProbeDir = mesh_.time().path()/faceProbeSubDir;
        }
        // Remove ".."
        faceProbeDir.clean();

        // ignore known fields, close streams for fields that no longer exist
        forAllIter(HashPtrTable<OFstream>, faceProbeFilePtrs_, iter)
        {
            if (!currentFields.erase(iter.key()))
            {
                if (debug)
                {
                    Info<< "close faceProbe stream: " << iter()->name() << endl;
                }

                delete faceProbeFilePtrs_.remove(iter);
            }
        }

        // currentFields now just has the new fields - open streams for them
        forAllConstIter(wordHashSet, currentFields, iter)
        {
            const word& fieldName = iter.key();

            // Create directory if does not exist.
            mkDir(faceProbeDir);

            OFstream* fPtr = new OFstream(faceProbeDir/fieldName);

            OFstream& fout = *fPtr;

            if (debug)
            {
                Info<< "open faceProbe stream: " << fout.name() << endl;
            }

            faceProbeFilePtrs_.insert(fieldName, fPtr);

            unsigned int w = IOstream::defaultPrecision() + 7;

            forAll(*this, faceProbei)
            {
                fout<< "# Probe " << faceProbei << ' ' << operator[](faceProbei)
                    << endl;
            }

            fout<< '#' << setw(IOstream::defaultPrecision() + 6)
                << "Probe";

            forAll(*this, faceProbei)
            {
                fout<< ' ' << setw(w) << faceProbei;
            }
            fout<< endl;

            fout<< '#' << setw(IOstream::defaultPrecision() + 6)
                << "Time" << endl;
        }
    }

    return nFields;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceProbes::faceProbes
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    pointField(0),
    mesh_
    (
        refCast<const fvMesh>
        (
            t.lookupObject<objectRegistry>
            (
                dict.lookupOrDefault("region", polyMesh::defaultRegion)
            )
        )
    ),
    loadFromFiles_(false),
    fieldSelection_(),
    fluxVector_(),
    fixedLocations_(true),
    interpolationScheme_("cell")
{
    read(dict);
}


Foam::faceProbes::faceProbes
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObject(name),
    pointField(0),
    mesh_(refCast<const fvMesh>(obr)),
    loadFromFiles_(loadFromFiles),
    fieldSelection_(),
    fluxVector_(),
    fixedLocations_(true),
    interpolationScheme_("cell")
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceProbes::~faceProbes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::faceProbes::read(const dictionary& dict)
{
    dict.lookup("faceProbeLocations") >> *this;
    dict.lookup("fields") >> fieldSelection_;
    dict.lookup("fluxVector") >> fluxVector_;

    dict.readIfPresent("fixedLocations", fixedLocations_);
    if (dict.readIfPresent("interpolationScheme", interpolationScheme_))
    {
        if (!fixedLocations_ && interpolationScheme_ != "cell")
        {
            WarningInFunction
                << "Only cell interpolation can be applied when "
                << "not using fixedLocations.  InterpolationScheme "
                << "entry will be ignored";
        }
    }

    // Initialise cells to sample from supplied locations
    findElements(mesh_);

    prepare();

    return true;
}


bool Foam::faceProbes::execute()
{
    return true;
}


bool Foam::faceProbes::write()
{
    if (size() && prepare())
    {
        sampleAndWrite(scalarFields_);
        sampleAndWrite(vectorFields_);
        sampleAndWrite(sphericalTensorFields_);
        sampleAndWrite(symmTensorFields_);
        sampleAndWrite(tensorFields_);

        sampleAndWriteSurfaceFields(surfaceScalarFields_);
        sampleAndWriteSurfaceFields(surfaceVectorFields_);
        sampleAndWriteSurfaceFields(surfaceSphericalTensorFields_);
        sampleAndWriteSurfaceFields(surfaceSymmTensorFields_);
        sampleAndWriteSurfaceFields(surfaceTensorFields_);
    }

    return true;
}


void Foam::faceProbes::updateMesh(const mapPolyMesh& mpm)
{
    DebugInfo<< "faceProbes: updateMesh" << endl;

    if (&mpm.mesh() != &mesh_)
    {
        return;
    }

    if (fixedLocations_)
    {
        findElements(mesh_);
    }
    else
    {
        if (debug)
        {
            Info<< "faceProbes: remapping sample locations" << endl;
        }

        // 1. Update cells
        {
            DynamicList<label> elems(elementList_.size());

            const labelList& reverseMap = mpm.reverseCellMap();
            forAll(elementList_, i)
            {
                label celli = elementList_[i];
                label newCelli = reverseMap[celli];
                if (newCelli == -1)
                {
                    // cell removed
                }
                else if (newCelli < -1)
                {
                    // cell merged
                    elems.append(-newCelli - 2);
                }
                else
                {
                    // valid new cell
                    elems.append(newCelli);
                }
            }

            elementList_.transfer(elems);
        }

        // 2. Update faces
        {
            DynamicList<label> elems(faceList_.size());

            const labelList& reverseMap = mpm.reverseFaceMap();
            forAll(faceList_, i)
            {
                label facei = faceList_[i];
                label newFacei = reverseMap[facei];
                if (newFacei == -1)
                {
                    // face removed
                }
                else if (newFacei < -1)
                {
                    // face merged
                    elems.append(-newFacei - 2);
                }
                else
                {
                    // valid new face
                    elems.append(newFacei);
                }
            }

            faceList_.transfer(elems);
        }
    }
}


void Foam::faceProbes::movePoints(const polyMesh& mesh)
{
    DebugInfo<< "faceProbes: movePoints" << endl;

    if (fixedLocations_ && &mesh == &mesh_)
    {
        findElements(mesh_);
    }
}


// ************************************************************************* //
