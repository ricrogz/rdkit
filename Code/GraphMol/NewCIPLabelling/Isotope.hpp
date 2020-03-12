//
//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once

#include <boost/unordered_map.hpp>
#include <vector>

namespace RDKit
{
namespace NewCIPLabelling
{

namespace
{

/**
 * From Blue Obelisk Data Repository.
 *
 * @see <a href="https://datahub.io/dataset/bodr">datahub.io</a>
 */
const auto isotopes = boost::unordered_map<std::pair<int, int>, double>({
    {{1, 0}, 1.007941},      // H
    {{1, 1}, 1.007825},      // H1
    {{1, 2}, 2.014102},      // H2
    {{2, 0}, 4.002602},      // He
    {{2, 3}, 3.016029},      // He3
    {{2, 4}, 4.002603},      // He4
    {{3, 0}, 6.940038},      // Li
    {{3, 6}, 6.015123},      // Li6
    {{3, 7}, 7.016005},      // Li7
    {{4, 0}, 9.012182},      // Be
    {{4, 9}, 9.012182},      // Be9
    {{5, 0}, 10.811028},     // B
    {{5, 10}, 10.012937},    // B10
    {{5, 11}, 11.009305},    // B11
    {{6, 0}, 12.010736},     // C
    {{6, 12}, 12.000000},    // C12
    {{6, 13}, 13.003355},    // C13
    {{7, 0}, 14.006703},     // N
    {{7, 14}, 14.003074},    // N14
    {{7, 15}, 15.000109},    // N15
    {{8, 0}, 15.999405},     // O
    {{8, 16}, 15.994915},    // O16
    {{8, 17}, 16.999132},    // O17
    {{8, 18}, 17.999161},    // O18
    {{9, 0}, 18.998403},     // F
    {{9, 19}, 18.998403},    // F19
    {{10, 0}, 20.180046},    // Ne
    {{10, 20}, 19.992440},   // Ne20
    {{10, 21}, 20.993847},   // Ne21
    {{10, 22}, 21.991385},   // Ne22
    {{11, 0}, 22.989769},    // Na
    {{11, 23}, 22.989769},   // Na23
    {{12, 0}, 24.305052},    // Mg
    {{12, 24}, 23.985042},   // Mg24
    {{12, 25}, 24.985837},   // Mg25
    {{12, 26}, 25.982593},   // Mg26
    {{13, 0}, 26.981539},    // Al
    {{13, 27}, 26.981539},   // Al27
    {{14, 0}, 28.085385},    // Si
    {{14, 28}, 27.976927},   // Si28
    {{14, 29}, 28.976495},   // Si29
    {{14, 30}, 29.973770},   // Si30
    {{15, 0}, 30.973762},    // P
    {{15, 31}, 30.973762},   // P31
    {{16, 0}, 32.066085},    // S
    {{16, 32}, 31.972071},   // S32
    {{16, 33}, 32.971459},   // S33
    {{16, 34}, 33.967867},   // S34
    {{16, 36}, 35.967081},   // S36
    {{17, 0}, 35.452938},    // Cl
    {{17, 35}, 34.968853},   // Cl35
    {{17, 37}, 36.965903},   // Cl37
    {{18, 0}, 39.947677},    // Ar
    {{18, 36}, 35.967545},   // Ar36
    {{18, 38}, 37.962732},   // Ar38
    {{18, 40}, 39.962383},   // Ar40
    {{19, 0}, 39.098301},    // K
    {{19, 39}, 38.963707},   // K39
    {{19, 40}, 39.963998},   // K40
    {{19, 41}, 40.961826},   // K41
    {{20, 0}, 40.078023},    // Ca
    {{20, 40}, 39.962591},   // Ca40
    {{20, 42}, 41.958618},   // Ca42
    {{20, 43}, 42.958767},   // Ca43
    {{20, 44}, 43.955482},   // Ca44
    {{20, 46}, 45.953693},   // Ca46
    {{20, 48}, 47.952534},   // Ca48
    {{21, 0}, 44.955912},    // Sc
    {{21, 45}, 44.955912},   // Sc45
    {{22, 0}, 47.866749},    // Ti
    {{22, 46}, 45.952632},   // Ti46
    {{22, 47}, 46.951763},   // Ti47
    {{22, 48}, 47.947946},   // Ti48
    {{22, 49}, 48.947870},   // Ti49
    {{22, 50}, 49.944791},   // Ti50
    {{23, 0}, 50.941467},    // V
    {{23, 50}, 49.947159},   // V50
    {{23, 51}, 50.943960},   // V51
    {{24, 0}, 51.996133},    // Cr
    {{24, 50}, 49.946044},   // Cr50
    {{24, 52}, 51.940508},   // Cr52
    {{24, 53}, 52.940649},   // Cr53
    {{24, 54}, 53.938880},   // Cr54
    {{25, 0}, 54.938045},    // Mn
    {{25, 55}, 54.938045},   // Mn55
    {{26, 0}, 55.845146},    // Fe
    {{26, 54}, 53.939611},   // Fe54
    {{26, 56}, 55.934938},   // Fe56
    {{26, 57}, 56.935394},   // Fe57
    {{26, 58}, 57.933276},   // Fe58
    {{27, 0}, 58.933195},    // Co
    {{27, 59}, 58.933195},   // Co59
    {{28, 0}, 58.693352},    // Ni
    {{28, 58}, 57.935343},   // Ni58
    {{28, 60}, 59.930786},   // Ni60
    {{28, 61}, 60.931056},   // Ni61
    {{28, 62}, 61.928345},   // Ni62
    {{28, 64}, 63.927966},   // Ni64
    {{29, 0}, 63.545640},    // Cu
    {{29, 63}, 62.929598},   // Cu63
    {{29, 65}, 64.927790},   // Cu65
    {{30, 0}, 65.395563},    // Zn
    {{30, 64}, 63.929142},   // Zn64
    {{30, 66}, 65.926033},   // Zn66
    {{30, 67}, 66.927127},   // Zn67
    {{30, 68}, 67.924844},   // Zn68
    {{30, 70}, 69.925319},   // Zn70
    {{31, 0}, 69.723066},    // Ga
    {{31, 69}, 68.925574},   // Ga69
    {{31, 71}, 70.924701},   // Ga71
    {{32, 0}, 72.612758},    // Ge
    {{32, 70}, 69.924247},   // Ge70
    {{32, 72}, 71.922076},   // Ge72
    {{32, 73}, 72.923459},   // Ge73
    {{32, 74}, 73.921178},   // Ge74
    {{32, 76}, 75.921403},   // Ge76
    {{33, 0}, 74.921597},    // As
    {{33, 75}, 74.921597},   // As75
    {{34, 0}, 78.959388},    // Se
    {{34, 74}, 73.922476},   // Se74
    {{34, 76}, 75.919214},   // Se76
    {{34, 77}, 76.919914},   // Se77
    {{34, 78}, 77.917309},   // Se78
    {{34, 80}, 79.916521},   // Se80
    {{34, 82}, 81.916699},   // Se82
    {{35, 0}, 79.903528},    // Br
    {{35, 79}, 78.918337},   // Br79
    {{35, 81}, 80.916291},   // Br81
    {{36, 0}, 83.799325},    // Kr
    {{36, 78}, 77.920365},   // Kr78
    {{36, 80}, 79.916379},   // Kr80
    {{36, 82}, 81.913484},   // Kr82
    {{36, 83}, 82.914136},   // Kr83
    {{36, 84}, 83.911507},   // Kr84
    {{36, 86}, 85.910611},   // Kr86
    {{37, 0}, 85.467664},    // Rb
    {{37, 85}, 84.911790},   // Rb85
    {{37, 87}, 86.909181},   // Rb87
    {{38, 0}, 87.616644},    // Sr
    {{38, 84}, 83.913425},   // Sr84
    {{38, 86}, 85.909260},   // Sr86
    {{38, 87}, 86.908877},   // Sr87
    {{38, 88}, 87.905612},   // Sr88
    {{39, 0}, 88.905848},    // Y
    {{39, 89}, 88.905848},   // Y89
    {{40, 0}, 91.223648},    // Zr
    {{40, 90}, 89.904704},   // Zr90
    {{40, 91}, 90.905646},   // Zr91
    {{40, 92}, 91.905041},   // Zr92
    {{40, 94}, 93.906315},   // Zr94
    {{40, 96}, 95.908273},   // Zr96
    {{41, 0}, 92.906378},    // Nb
    {{41, 93}, 92.906378},   // Nb93
    {{42, 0}, 95.931292},    // Mo
    {{42, 92}, 91.906811},   // Mo92
    {{42, 94}, 93.905088},   // Mo94
    {{42, 95}, 94.905842},   // Mo95
    {{42, 96}, 95.904680},   // Mo96
    {{42, 97}, 96.906022},   // Mo97
    {{42, 98}, 97.905408},   // Mo98
    {{42, 100}, 99.907477},  // Mo100
    {{44, 0}, 101.064945},   // Ru
    {{44, 96}, 95.907598},   // Ru96
    {{44, 98}, 97.905287},   // Ru98
    {{44, 99}, 98.905939},   // Ru99
    {{44, 100}, 99.904220},  // Ru100
    {{44, 101}, 100.905582}, // Ru101
    {{44, 102}, 101.904349}, // Ru102
    {{44, 104}, 103.905433}, // Ru104
    {{45, 0}, 102.905504},   // Rh
    {{45, 103}, 102.905504}, // Rh103
    {{46, 0}, 106.415329},   // Pd
    {{46, 102}, 101.905609}, // Pd102
    {{46, 104}, 103.904036}, // Pd104
    {{46, 105}, 104.905085}, // Pd105
    {{46, 106}, 105.903486}, // Pd106
    {{46, 108}, 107.903892}, // Pd108
    {{46, 110}, 109.905153}, // Pd110
    {{47, 0}, 107.868151},   // Ag
    {{47, 107}, 106.905097}, // Ag107
    {{47, 109}, 108.904752}, // Ag109
    {{48, 0}, 112.411552},   // Cd
    {{48, 106}, 105.906459}, // Cd106
    {{48, 108}, 107.904184}, // Cd108
    {{48, 110}, 109.903002}, // Cd110
    {{48, 111}, 110.904178}, // Cd111
    {{48, 112}, 111.902758}, // Cd112
    {{48, 113}, 112.904402}, // Cd113
    {{48, 114}, 113.903359}, // Cd114
    {{48, 116}, 115.904756}, // Cd116
    {{49, 0}, 114.818086},   // In
    {{49, 113}, 112.904058}, // In113
    {{49, 115}, 114.903878}, // In115
    {{50, 0}, 118.710108},   // Sn
    {{50, 112}, 111.904818}, // Sn112
    {{50, 114}, 113.902779}, // Sn114
    {{50, 115}, 114.903342}, // Sn115
    {{50, 116}, 115.901741}, // Sn116
    {{50, 117}, 116.902952}, // Sn117
    {{50, 118}, 117.901603}, // Sn118
    {{50, 119}, 118.903308}, // Sn119
    {{50, 120}, 119.902195}, // Sn120
    {{50, 122}, 121.903439}, // Sn122
    {{50, 124}, 123.905274}, // Sn124
    {{51, 0}, 121.759786},   // Sb
    {{51, 121}, 120.903816}, // Sb121
    {{51, 123}, 122.904214}, // Sb123
    {{52, 0}, 127.603128},   // Te
    {{52, 120}, 119.904020}, // Te120
    {{52, 122}, 121.903044}, // Te122
    {{52, 123}, 122.904270}, // Te123
    {{52, 124}, 123.902818}, // Te124
    {{52, 125}, 124.904431}, // Te125
    {{52, 126}, 125.903312}, // Te126
    {{52, 128}, 127.904463}, // Te128
    {{52, 130}, 129.906224}, // Te130
    {{53, 0}, 126.904473},   // I
    {{53, 127}, 126.904473}, // I127
    {{54, 0}, 131.292481},   // Xe
    {{54, 124}, 123.905893}, // Xe124
    {{54, 126}, 125.904274}, // Xe126
    {{54, 128}, 127.903531}, // Xe128
    {{54, 129}, 128.904779}, // Xe129
    {{54, 130}, 129.903508}, // Xe130
    {{54, 131}, 130.905082}, // Xe131
    {{54, 132}, 131.904154}, // Xe132
    {{54, 134}, 133.905395}, // Xe134
    {{54, 136}, 135.907219}, // Xe136
    {{55, 0}, 132.905452},   // Cs
    {{55, 133}, 132.905452}, // Cs133
    {{56, 0}, 137.326892},   // Ba
    {{56, 130}, 129.906321}, // Ba130
    {{56, 132}, 131.905061}, // Ba132
    {{56, 134}, 133.904508}, // Ba134
    {{56, 135}, 134.905689}, // Ba135
    {{56, 136}, 135.904576}, // Ba136
    {{56, 137}, 136.905827}, // Ba137
    {{56, 138}, 137.905247}, // Ba138
    {{57, 0}, 138.905454},   // La
    {{57, 138}, 137.907112}, // La138
    {{57, 139}, 138.906353}, // La139
    {{58, 0}, 140.115726},   // Ce
    {{58, 136}, 135.907172}, // Ce136
    {{58, 138}, 137.905991}, // Ce138
    {{58, 140}, 139.905439}, // Ce140
    {{58, 142}, 141.909244}, // Ce142
    {{59, 0}, 140.907653},   // Pr
    {{59, 141}, 140.907653}, // Pr141
    {{60, 0}, 144.236131},   // Nd
    {{60, 142}, 141.907723}, // Nd142
    {{60, 143}, 142.909814}, // Nd143
    {{60, 144}, 143.910087}, // Nd144
    {{60, 145}, 144.912574}, // Nd145
    {{60, 146}, 145.913117}, // Nd146
    {{60, 148}, 147.916893}, // Nd148
    {{60, 150}, 149.920891}, // Nd150
    {{62, 0}, 150.366349},   // Sm
    {{62, 144}, 143.911999}, // Sm144
    {{62, 147}, 146.914898}, // Sm147
    {{62, 148}, 147.914823}, // Sm148
    {{62, 149}, 148.917185}, // Sm149
    {{62, 150}, 149.917276}, // Sm150
    {{62, 152}, 151.919732}, // Sm152
    {{62, 154}, 153.922209}, // Sm154
    {{63, 0}, 151.964370},   // Eu
    {{63, 151}, 150.919850}, // Eu151
    {{63, 153}, 152.921230}, // Eu153
    {{64, 0}, 157.252122},   // Gd
    {{64, 152}, 151.919791}, // Gd152
    {{64, 154}, 153.920866}, // Gd154
    {{64, 155}, 154.922622}, // Gd155
    {{64, 156}, 155.922123}, // Gd156
    {{64, 157}, 156.923960}, // Gd157
    {{64, 158}, 157.924104}, // Gd158
    {{64, 160}, 159.927054}, // Gd160
    {{65, 0}, 158.925347},   // Tb
    {{65, 159}, 158.925347}, // Tb159
    {{66, 0}, 162.497034},   // Dy
    {{66, 156}, 155.924283}, // Dy156
    {{66, 158}, 157.924409}, // Dy158
    {{66, 160}, 159.925198}, // Dy160
    {{66, 161}, 160.926933}, // Dy161
    {{66, 162}, 161.926798}, // Dy162
    {{66, 163}, 162.928731}, // Dy163
    {{66, 164}, 163.929175}, // Dy164
    {{67, 0}, 164.930322},   // Ho
    {{67, 165}, 164.930322}, // Ho165
    {{68, 0}, 167.256304},   // Er
    {{68, 162}, 161.928778}, // Er162
    {{68, 164}, 163.929200}, // Er164
    {{68, 166}, 165.930293}, // Er166
    {{68, 167}, 166.932048}, // Er167
    {{68, 168}, 167.932370}, // Er168
    {{68, 170}, 169.935464}, // Er170
    {{69, 0}, 168.934213},   // Tm
    {{69, 169}, 168.934213}, // Tm169
    {{70, 0}, 173.037696},   // Yb
    {{70, 168}, 167.933897}, // Yb168
    {{70, 170}, 169.934762}, // Yb170
    {{70, 171}, 170.936326}, // Yb171
    {{70, 172}, 171.936382}, // Yb172
    {{70, 173}, 172.938211}, // Yb173
    {{70, 174}, 173.938862}, // Yb174
    {{70, 176}, 175.942572}, // Yb176
    {{71, 0}, 174.966721},   // Lu
    {{71, 175}, 174.940772}, // Lu175
    {{71, 176}, 175.942686}, // Lu176
    {{72, 0}, 178.484972},   // Hf
    {{72, 174}, 173.940046}, // Hf174
    {{72, 176}, 175.941409}, // Hf176
    {{72, 177}, 176.943221}, // Hf177
    {{72, 178}, 177.943699}, // Hf178
    {{72, 179}, 178.945816}, // Hf179
    {{72, 180}, 179.946550}, // Hf180
    {{73, 0}, 180.947876},   // Ta
    {{73, 180}, 179.947465}, // Ta180
    {{73, 181}, 180.947996}, // Ta181
    {{74, 0}, 183.841778},   // W
    {{74, 180}, 179.946704}, // W180
    {{74, 182}, 181.948204}, // W182
    {{74, 183}, 182.950223}, // W183
    {{74, 184}, 183.950931}, // W184
    {{74, 186}, 185.954364}, // W186
    {{75, 0}, 186.206707},   // Re
    {{75, 185}, 184.952955}, // Re185
    {{75, 187}, 186.955753}, // Re187
    {{76, 0}, 190.224863},   // Os
    {{76, 184}, 183.952489}, // Os184
    {{76, 186}, 185.953838}, // Os186
    {{76, 187}, 186.955751}, // Os187
    {{76, 188}, 187.955838}, // Os188
    {{76, 189}, 188.958148}, // Os189
    {{76, 190}, 189.958447}, // Os190
    {{76, 192}, 191.961481}, // Os192
    {{77, 0}, 192.216056},   // Ir
    {{77, 191}, 190.960594}, // Ir191
    {{77, 193}, 192.962926}, // Ir193
    {{78, 0}, 195.077808},   // Pt
    {{78, 190}, 189.959932}, // Pt190
    {{78, 192}, 191.961038}, // Pt192
    {{78, 194}, 193.962680}, // Pt194
    {{78, 195}, 194.964791}, // Pt195
    {{78, 196}, 195.964952}, // Pt196
    {{78, 198}, 197.967893}, // Pt198
    {{79, 0}, 196.966569},   // Au
    {{79, 197}, 196.966569}, // Au197
    {{80, 0}, 200.599167},   // Hg
    {{80, 196}, 195.965833}, // Hg196
    {{80, 198}, 197.966769}, // Hg198
    {{80, 199}, 198.968280}, // Hg199
    {{80, 200}, 199.968326}, // Hg200
    {{80, 201}, 200.970302}, // Hg201
    {{80, 202}, 201.970643}, // Hg202
    {{80, 204}, 203.973494}, // Hg204
    {{81, 0}, 204.383332},   // Tl
    {{81, 203}, 202.972344}, // Tl203
    {{81, 205}, 204.974428}, // Tl205
    {{82, 0}, 207.216908},   // Pb
    {{82, 204}, 203.973044}, // Pb204
    {{82, 206}, 205.974465}, // Pb206
    {{82, 207}, 206.975897}, // Pb207
    {{82, 208}, 207.976652}, // Pb208
    {{83, 0}, 208.980399},   // Bi
    {{83, 209}, 208.980399}, // Bi209
    {{90, 0}, 232.038055},   // Th
    {{90, 232}, 232.038055}, // Th232
    {{91, 0}, 231.035884},   // Pa
    {{91, 231}, 231.035884}, // Pa231
    {{92, 0}, 238.028918},   // U
    {{92, 234}, 234.040952}, // U234
    {{92, 235}, 235.043930}, // U235
    {{92, 238}, 238.050788}, // U238
});
} // namespace

namespace Isotope
{
static double findAtomicWeight(int atomnum, int massnum)
{
    try {
        return isotopes.at({atomnum, massnum});
    } catch (const std::out_of_range&) {
        return static_cast<double>(massnum);
    }
}
} // namespace Isotope

} // namespace NewCIPLabelling
} // namespace RDKit
