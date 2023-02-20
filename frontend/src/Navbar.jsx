import React from "react";
import { Container, Row, Col, Form, Navbar, Nav } from "react-bootstrap";
export default function Navigation() {
    return (
        <Navbar bg="dark" expand="lg" variant="dark">
            <Container>
                <Navbar.Brand href="/">Lithium</Navbar.Brand>
                <Navbar.Toggle aria-controls="basic-navbar-nav" />
                <Navbar.Collapse id="basic-navbar-nav">
                    <Nav className="me-auto">
                        <Nav.Link href="/substances">Substances</Nav.Link>
                        <Nav.Link href="/similarity">Similarity</Nav.Link>
                        <Nav.Link href="/substructure">Substructure</Nav.Link>
                    </Nav>
                </Navbar.Collapse>
            </Container>
        </Navbar>
    );
}