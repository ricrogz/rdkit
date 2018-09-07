# SGroups in RDKit
Ricardo Rodriguez-Schmidt

September 2018

*This document is a DRAFT*

## Overview

MDL SD files suport SGroups as a way to store some types of extendend intra-molecular information. Some usage examples of these SGroup are the following:

1. Identification of repeating groups in polymers, i.e. identification of monomer units and, in case of copolymers, how these monomers distribute in the polymer chain.
1. Labeling of relevant sections of a molecule.
1. Detailing depiction information for the molecule or other defined SGroups.
1. Storing information relative to parts of the molecule in the form of data fields.

This document describes an intent to make a `SGroup` data structure available to RDKit molecules with the goal of allowing parsing, storing, serializing and writing back this information, with the goal to achieve  conversion of MDL SD files that contain SGroups using RDkit without loss of information. At this point, no further intents will be made to allow or support manipulation or depiction of the information contained in SGroups, but this might be added in the future.

## Documentation

As a reference for the design and implementation of the `SGroup` data class, the following document has been used:

[http://infochim.u-strasbg.fr/recherche/Download/Fragmentor/MDL_SDF.pdf](http://infochim.u-strasbg.fr/recherche/Download/Fragmentor/MDL_SDF.pdf)

The document describes the syntax for specifying SGRoup information both in version V2000 and V3000 of MDL SD files. Both of these have been implemented for input and output, so that interconversion is possible. *Note*: All SGroup available to V2000 and V3000 have been implemented with the exception of the V2000 "CRS" / V3000 "XBHEAD/XBCORR" labels.

Unfortunately, the available documentation is somewhat incomplete at some points, and does not provide global a overview of the SGroup syntax, resultin in some details being specified only for V2000 or V3000 format, and not for the other one. Best effort has been made for compatibility.

## Testing

Due to the poor availability of molecules with SGroup information, the implementation has only been lightly tested with some samples already included with RDKit that have simple SGroup information, as well as with custom tailored molecules (which do not make real chemical sense).

This means that the most of the deductions and assumptions made from the documentation have not been tested.

## Proposed SGroup Data Structure

The SGroup data structure is implemented through the self contained `SGroup` class. `SGroup` objects are not meant to exist independently, but only as objects owned by `ROMol` (and derivate) objects. In these, they are stored inside a vector of shared pointers (`std::vector<boost::shared_ptr<SGroup>>`).

The proposed implementation of the SGroup class has the following private properties, together with getter/setter methods associated to them. References to SGroup format specification are relative to the V3000 version. Most of these may be unset, empty or have default values if not specified.

`unsigned int` **d_id**: Storage for the *extindex* value. The *index* value is not explicitly stored, as it is meant for reference inside files, and is made redundant by the containing vector's indexing. May be set or have the value 0 (meant to represent automatic assignment).

`unsigned int` **d_compno**: Storage for the *COMPNO* value. Default is 0, which is ignored.

`string` **d_type**: Mandatory. Type of the *SGroup*, represented by a three letter string. This is validated on input, and some of the behavior of the SGroup depends on the declared type.

`ROMol \*` **dp_mol**: Pointer to the owning `ROMol` object. Should contain a `nullptr` only at creation of the object.

`SGroup \*` **d_parent**: Pointer to the parent `SGroup` object. `nullptr` if the property is not set.

`vector<Atom \*>` **d_atoms**: Atoms contained in the `SGroup`.

`vector<Atom \*>` **d_patoms**: PATOMS contained in the `SGroup`.

`vector<Bond \*>` **d_bonds**: Bonds in the SGroup. At least one of the atoms at its ends should belong to the SGroup (not enforced).

`vector<Bracket>` **d_brackets**: Storage for BRKXYZ. Documentation does not specify how many brackets a `SGroup` may contain, although it seems reasonable to have 2.

`vector<CState>` **d_cstates**: Storage for CSTATEs. Documentation does not specify how many CSTATEs a `SGroup` may contain.

`vector<string>` **d_dataFields**: Data fields associated to the SGroup. There may be multiple fields for the same `SGroup`.

`vector<AttachPoint>` **d_saps**: Storage for SAPs. Documentation does not specify how many SAPs a `SGroup` may contain.

`unordered_map<string, string>` **d_strProp**: Storage for all other SGroup labels that are stored as text strings. Non-standard labels in V3000 format may also be stored in this container, if required.

Besides these properties, the `SGroup` class also defines three data types:

* **Bracket**:
    ```
    typedef array<RDGeom::Point3D, 3> Bracket
    ```

Coordinates for the Brackets to be drawn around the SGroup. V2000 only allows for 2 points with *x* and *y* coordinates, although V3000 also has *z* and a third point. For compatibility between formats, *z* and the third point are all zeros.

* **AttachPoint**:
    ```
    struct AttachPoint {
        Atom *aAtom
        Atom *lvAtom
        string id
    };
    ```

Indicates the *aAtom* to which the atom where the SGroup is supposed to be attached. *lvAtom* indicates an atom leaving the molecule when the SGroup attaches. If no atom is leaving, lvAtom will be a `nullptr`. It also can be equal to *aAtom* (in V2000, having the same index; in V3000 by indicating 'aidx' in the file).

* **CState**:
    ```
    struct CState {
        Bond *bond;
        shared_ptr<RDGeom::Point3D> vector;
    };
    ```

Depending on the type of the SGroup, *CState* may only contain a *Bond* (general behavior; vector is `nullptr`), or also a vector, in case of abbreviation SGroups (type *SUP*? -- deduced from documentation, may be wrong).

**Remark**:

The `SGroup` class contains two private "remapping" methods:

```
//! Update atoms, patoms, bonds with the ones at matching indexes in other_mol
void remap_atoms_bonds_to_new_mol(ROMol *other_mol);

//! Update parent sgroup with the one at the same index in other_mol
void remap_parent_sgroup_to_new_mol(ROMol *other_mol);
```

These are meant to be called after copying and/or assigning an `SGroup` to a new molecule.

The first one updates the pointers to the Atoms and Bonds contained in the `SGroup` with those at the same indexes in the new owning molecule. It can be called at creation of the `Sgroup` since the Atoms and Bonds should already be present on the new owning molecule.

The second method updates the pointer to the parent `SGroup`, and should be called only once all `SGroups` have been created/copied on the new molecule, since a `SGroup` can be child to another one that is located after it in the molecule's `d_sgroups` vector (and therefore, might not be initialized yet).

