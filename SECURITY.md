# Security Policy

## Supported Versions

We provide security updates for the following versions of MutationScan:

| Version | Supported          |
| ------- | ------------------ |
| 1.0.x   | :white_check_mark: |
| < 1.0   | :x:                |

## Reporting a Vulnerability

The MutationScan team takes security vulnerabilities seriously. We appreciate your efforts to responsibly disclose any security concerns.

### How to Report

**Please do not report security vulnerabilities through public GitHub issues.**

Instead, please send a detailed report to: **vihaankulkarni29@gmail.com**

Include the following information:
- Description of the vulnerability
- Steps to reproduce the issue
- Potential impact assessment
- Suggested fix (if you have one)

### What to Expect

1. **Acknowledgment**: We will acknowledge receipt of your vulnerability report within 48 hours
2. **Assessment**: We will assess the vulnerability and determine its severity within 7 days
3. **Updates**: We will keep you informed of our progress toward resolving the issue
4. **Resolution**: We aim to resolve security vulnerabilities within 30 days of confirmation

### Security Considerations

MutationScan processes genomic data and interacts with external APIs. Key security considerations include:

#### Data Handling
- **Input Validation**: All user inputs are validated and sanitized
- **File Processing**: Genomic files are processed in sandboxed environments
- **Temporary Files**: Temporary files are securely deleted after processing

#### External Communications
- **NCBI API**: All API communications use HTTPS
- **Rate Limiting**: Built-in rate limiting prevents abuse
- **Authentication**: No sensitive credentials are stored in code

#### Dependencies
- **Regular Updates**: Dependencies are regularly updated for security patches
- **Vulnerability Scanning**: Automated scanning for known vulnerabilities
- **Minimal Dependencies**: We minimize third-party dependencies to reduce attack surface

### Scope

This security policy applies to:
- Core MutationScan pipeline components
- Command-line interfaces
- API integrations
- Documentation and example code

### Out of Scope

The following are generally not considered security vulnerabilities:
- Issues requiring physical access to the host system
- Social engineering attacks
- Vulnerabilities in third-party dependencies (please report to upstream)
- Denial of service through excessive computational load (expected behavior)

### Recognition

We believe in recognizing security researchers who help improve our software. With your permission, we will:
- Credit you in our security advisories
- List you in our acknowledgments
- Provide a reference for responsible disclosure

Thank you for helping keep MutationScan and our users safe!