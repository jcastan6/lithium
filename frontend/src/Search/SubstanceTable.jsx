import React from "react";
import { useEffect, useImperativeHandle, forwardRef } from "react";
import { Col, Row, Button, Card, Table, Pagination, Navbar, Spinner, Dropdown } from 'react-bootstrap';
import MoleculeStructure from '../RDKit/MoleculeStructure';
import CsvDownloadButton from 'react-json-to-csv'

const SubstanceTable = forwardRef((props, ref) => {
    const [molecules, setMolecules] = React.useState([]);
    const [substructure, setSubstructure] = React.useState(props.subStructure);

    const [page, setPage] = React.useState(1);
    const [perPage, setPerPage] = React.useState(50);
    const [loading, setLoading] = React.useState(false);
    const [total, setTotal] = React.useState(0);

    useImperativeHandle(ref, () => ({
        getMolecules(substructure = null) {
            getMolecules(substructure);
        }
    }));

    useEffect(() => {
        getTotal();
    }, []);

    useEffect(() => {
        getMolecules();
    }, [page, perPage]);


    function getMolecules(substructure = null) {

        setLoading(true);
        fetch(props.url + "?page=" + page + "&per_page=" + perPage)
            .then(response => response.json())
            .then(data => {

                setMolecules(data)
                setLoading(false);
            });
    }

    function getTotal() {
        fetch(props.url + "/total")
            .then(response => response.json())
            .then(data => {
                setTotal(data)
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
        if (page > Math.ceil(total / perPage) - 2) {
            start = Math.ceil(total / perPage) - 4;
            end = Math.ceil(total / perPage);
        }
        //if total pages is less than 5, show all pages
        if (Math.ceil(total / perPage) < 5) {
            start = 1;
            end = Math.ceil(total / perPage);
        }
        //if start is less than 1, set to 1
        if (start < 1) {
            start = 1;
        }
        //if end is greater than total pages, set to total pages
        if (end > Math.ceil(total / perPage)) {
            end = Math.ceil(total / perPage);
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
                        <Pagination.Next onClick={() => page >= Math.ceil(total / perPage) ? setPage(Math.ceil(total / perPage)) : setPage(page + 1)} />
                        <Pagination.Last onClick={() => setPage(Math.ceil(total / perPage))} />

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
                        <Button>{total} results</Button>
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
                <Row>
                    {loading && <Spinner animation="border" role="status" className="mx-auto">
                        <span className="visually-hidden">Loading...</span> </Spinner>}
                    {!loading && total > 0 &&
                        // display molecules in cards
                        molecules.map((molecule, index) => {

                            return (
                                <Col lg={2}>
                                    <Card key={index} className="mb-3">
                                        <Card.Body>
                                            <Card.Title className="text-center">{molecule.id}</Card.Title>
                                            <MoleculeStructure
                                                structure={molecule.smiles}
                                                svgMode
                                            />

                                        </Card.Body>
                                    </Card>
                                </Col>
                            )
                        })
                    }

                </Row>
            </Card.Body>
        </Card >
    )

})
export default SubstanceTable;
