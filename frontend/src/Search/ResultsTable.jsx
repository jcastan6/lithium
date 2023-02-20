import React from "react";
import { useEffect, useImperativeHandle, forwardRef } from "react";
import { Col, Button, Card, Table, Pagination, Navbar, Spinner, Dropdown } from 'react-bootstrap';
import MoleculeStructure from '../RDKit/MoleculeStructure';
import CsvDownloadButton from 'react-json-to-csv'

const ResultsTable = forwardRef((props, ref) => {
    const [molecules, setMolecules] = React.useState([]);
    const [substructure, setSubstructure] = React.useState(props.subStructure);

    const [page, setPage] = React.useState(1);
    const [perPage, setPerPage] = React.useState(20);
    const [loading, setLoading] = React.useState(false);
    const [total, setTotal] = React.useState(0);
    const [totalPages, setTotalPages] = React.useState(0);

    useImperativeHandle(ref, () => ({
        getMolecules(substructure = null) {
            getMolecules(substructure);
        }
    }));





    function getMolecules(substructure = null) {

        setLoading(true);
        fetch(props.url + "?smiles=" + props.smiles)
            .then(response => response.json())
            .then(data => {
                if (substructure) {
                    setSubstructure(substructure);
                }

                setMolecules(data)
                setLoading(false);
            });

    }

    function buildPagination() {
        let pagination = [];
        //show 2 pages before and after current page
        let start = page - 2;
        let end = page + 2;
        //if current page is less than 3, show first 5 pages
        if (page < 3) {
            start = 1;
            end = 5;
        }
        //if current page is greater than total pages - 2, show last 5 pages
        if (page > Math.ceil(molecules.length / perPage) - 2) {
            start = Math.ceil(molecules.length / perPage) - 4;
            end = Math.ceil(molecules.length / perPage);
        }
        //if total pages is less than 5, show all pages
        if (Math.ceil(molecules.length / perPage) < 5) {
            start = 1;
            end = Math.ceil(molecules.length / perPage);
        }
        //if start is less than 1, set to 1
        if (start < 1) {
            start = 1;
        }
        //if end is greater than total pages, set to total pages
        if (end > Math.ceil(molecules.length / perPage)) {
            end = Math.ceil(molecules.length / perPage);
        }
        //build pagination
        for (let i = start; i <= end; i++) {
            pagination.push(
                <Pagination.Item key={i} active={i === page} onClick={() => setPage(i)}>
                    {i}
                </Pagination.Item>,
            );
        }


        return pagination;
    }

    useEffect(() => {
        setPage(1);
    }, [perPage]);

    return (
        <Card>
            <Card.Header>Results</Card.Header>
            <Card.Body>
                <Navbar bg="clear" className=''>
                    <Pagination style={{ "marginBottom": "auto", }}>
                        <Pagination.First onClick={() => setPage(1)} />
                        <Pagination.Prev onClick={() => page <= 1 ? setPage(1) : setPage(page - 1)} />
                        {buildPagination()}
                        <Pagination.Next onClick={() => page >= Math.ceil(molecules.length / perPage) ? setPage(Math.ceil(molecules.length / perPage)) : setPage(page + 1)} />
                        <Pagination.Last onClick={() => setPage(Math.ceil(molecules.length / perPage))} />

                    </Pagination>
                    &nbsp;
                    <Dropdown>
                        <Dropdown.Toggle variant="" id="dropdown-basic">
                            {perPage}
                        </Dropdown.Toggle>
                        <Dropdown.Menu>
                            <Dropdown.Item onClick={() => setPerPage(10)}>10</Dropdown.Item>
                            <Dropdown.Item onClick={() => setPerPage(20)}>20</Dropdown.Item>
                            <Dropdown.Item onClick={() => setPerPage(50)}>50</Dropdown.Item>
                            <Dropdown.Item onClick={() => setPerPage(100)}>100</Dropdown.Item>
                        </Dropdown.Menu>
                    </Dropdown>
                    &nbsp;



                    <Navbar.Collapse className='justify-content-end align-middle'>
                        <Button>{molecules.length} results</Button>
                        &nbsp;
                        <Dropdown>
                            <Dropdown.Toggle variant="success" id="dropdown-basic">
                                Download
                            </Dropdown.Toggle>

                            <Dropdown.Menu>
                                <CsvDownloadButton data={molecules} filename='search_results.csv' className='dropdown-item' target='_blank' delimiter=",">CSV</CsvDownloadButton>
                                <Dropdown.Item href={`data:text/json;charset=utf-8,${encodeURIComponent(
                                    JSON.stringify(molecules)
                                )}`} download='search_results.json'>JSON</Dropdown.Item>

                            </Dropdown.Menu>


                        </Dropdown>

                    </Navbar.Collapse>
                </Navbar>

                {loading && <Spinner animation="border" role="status" className="mx-auto"><span className="visually-hidden">Loading...</span> </Spinner>}
                {!loading && molecules.length > 0 &&
                    <Table striped bordered hover>
                        <thead>
                            <tr>
                                <th>Structure</th>

                                {molecules.length > 0 && Object.keys(molecules[0]).map((key, index) => (
                                    <th key={index}>{key}</th>
                                ))}
                            </tr>
                        </thead>
                        <tbody>
                            {molecules.slice((page - 1) * perPage, page * perPage).map((molecule, index) => (
                                <tr key={index}>
                                    <td><MoleculeStructure id={molecule.smiles} structure={molecule.smiles} height={150} width={150} svgMode
                                        subStructure={substructure}
                                    /></td>
                                    {Object.keys(molecule).map((key, index) => (
                                        <td key={index}>{molecule[key]}</td>
                                    ))}
                                </tr>
                            ))}
                        </tbody>
                    </Table>
                }
            </Card.Body>
        </Card>
    )

})
export default ResultsTable;
