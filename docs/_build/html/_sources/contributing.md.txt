# Contributing

Contributions to radiosonde-teamxuk are welcome! This guide will help you get started.

## Getting Started

### Development Setup

1. Fork the repository on GitHub
2. Clone your fork locally:
   ```bash
   git clone https://github.com/YOUR_USERNAME/radiosonde-teamxuk.git
   cd radiosonde-teamxuk
   ```

3. Install in editable mode with development dependencies:
   ```bash
   pip install -e .
   ```

4. Create a new branch for your feature or bugfix:
   ```bash
   git checkout -b feature/your-feature-name
   ```

## Development Guidelines

### Code Style

- Follow PEP 8 Python style guidelines
- Use meaningful variable and function names
- Add docstrings to all functions and classes
- Keep functions focused and modular

### Documentation

- Update documentation for any new features
- Include examples in docstrings
- Update the changelog for significant changes

### Testing

Before submitting changes:

1. Test your changes with sample data:
   ```bash
   # Test processing
   process-teamxuk-radiosondes test_edt_output/ test_output/
   
   # Test quicklooks
   generate-teamxuk-quicklooks test_output/ --stability
   ```

2. Verify the package version displays correctly:
   ```bash
   process-teamxuk-radiosondes --version
   ```

3. Check that imports work correctly:
   ```python
   from radiosonde_teamxuk import process_radiosondes
   from radiosonde_teamxuk import generate_quicklooks
   ```

## Submitting Changes

### Pull Request Process

1. Commit your changes with clear commit messages:
   ```bash
   git add .
   git commit -m "Add feature: description of your changes"
   ```

2. Push to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```

3. Open a Pull Request on GitHub:
   - Provide a clear title and description
   - Reference any related issues
   - Include screenshots for visual changes
   - List any breaking changes

### Pull Request Checklist

- [ ] Code follows the project style
- [ ] Documentation has been updated
- [ ] Changes have been tested with sample data
- [ ] Commit messages are clear and descriptive
- [ ] No unnecessary files are included

## Types of Contributions

### Bug Reports

If you find a bug:

1. Check if it's already reported in Issues
2. Create a new issue with:
   - Clear description of the problem
   - Steps to reproduce
   - Expected vs. actual behavior
   - Your environment (Python version, OS, etc.)
   - Sample data or code if possible

### Feature Requests

For new features:

1. Open an issue to discuss the feature first
2. Explain the use case and benefits
3. Consider implementation approach
4. Wait for feedback before starting work

### Code Contributions

Areas where contributions are welcome:

- **Additional derived variables**: Add more meteorological calculations
- **Plotting enhancements**: Improve or add new visualization options
- **Data format support**: Add support for other radiosonde formats
- **Performance improvements**: Optimize data processing
- **Error handling**: Improve robustness and error messages
- **Testing**: Add unit tests and integration tests
- **Documentation**: Improve examples and tutorials

### Documentation Improvements

Documentation can always be improved:

- Fix typos or unclear sections
- Add more examples
- Improve API documentation
- Create tutorials for specific use cases

## Code of Conduct

### Our Standards

- Be respectful and inclusive
- Welcome newcomers
- Focus on constructive feedback
- Assume good intentions

### Unacceptable Behavior

- Harassment or discriminatory language
- Trolling or insulting comments
- Publishing others' private information
- Other unprofessional conduct

## Questions?

If you have questions about contributing:

- Open an issue for discussion
- Check existing issues and pull requests
- Review the documentation

## License

By contributing, you agree that your contributions will be licensed under the same license as the project (see LICENSE file).

## Recognition

Contributors will be acknowledged in:
- The project README
- Release notes for significant contributions
- Git commit history

Thank you for contributing to radiosonde-teamxuk!
